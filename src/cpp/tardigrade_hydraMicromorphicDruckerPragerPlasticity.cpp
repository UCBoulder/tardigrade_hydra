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

        variableType weakMac( const variableType &x, const variableType &a ){
            /*!
             * Compute a weakened Macaulay bracket based on the sigmoid function
             * 
             * \param &x: The value to compute the weakened sigmoid function of
             * \param &a: The a parameter that determines how weak the bracket is
             */

            return std::log( std::exp( a * x ) + 1 ) / a;

        }

        variableType weakMac( const variableType &x, const variableType &a, variableType &dmacdx ){
            /*!
             * Compute a weakened Macaulay bracket based on the sigmoid function
             * 
             * \param &x: The value to compute the weakened sigmoid function of
             * \param &a: The a parameter that determines how weak the bracket is
             * \param &dmacdx: The derivative of the bracket w.r.t. x
             */

            dmacdx = 1. / ( 1. + std::exp( -a * x ) );

            return std::log( std::exp( a * x ) + 1 ) / a;

        }

        void computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B ){
            /*!
             * Compute the Drucker-Prager internal parameters
             *
             * \param &frictionAngle: The material friction angle ( 0 < frictionAngle < pi / 2 );
             * \param &beta: The beta parameter.
             * \param &A: The A parameter.
             * \param &B: The B parameter.
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
             * \param &stressMeasure: The stress measure
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The deformation gradients preceding the configuration of the stress measure
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
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
             * \param &referenceStressMeasure: The stress measure in the reference configuration
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The preceding deformation gradient.
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
             * \param &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
             * \param &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
             * \param &dFdPrecedingF: The Jacobian of the yield surface w.r.t. the preceding deformation gradient from the stress-measure's configuration
             * \param tol: The tolerance used to prevent nans in the Jacobians
             */

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );
   
            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatVector dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );
 
            //Compute the decomposition of the stress
            variableType pressure;
            variableVector deviatoricReferenceStress;
    
            variableVector dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );
   
            variableVector dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim );
            variableVector dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG, dRCGdPrecedingF, 1, sot_dim, sot_dim, sot_dim );
 
            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );
    
            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Evaluate the jacobians
            variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );
    
            dFdStress = tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdStress, 1, sot_dim, sot_dim, sot_dim )
                      + BAngle * dPressuredStress;
    
            dFdc = - AAngle;
    
            dFdPrecedingF = tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdPrecedingF, 1, sot_dim, sot_dim, sot_dim )
                          + BAngle * dPressuredPrecedingF;
        }

        void computeSecondOrderDruckerPragerYieldEquation( const variableVector &stressMeasure, const variableType &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdPrecedingF, variableVector &d2FdStress2,
                                                           variableVector &d2FdStressdPrecedingF, double tol ){
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
             * \param &stressMeasure: The stress measure
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The preceding deformation gradient
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
             * \param &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
             * \param &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
             * \param &dFdPrecedingF: The Jacobian of the yield surface w.r.t. the preceding
             *     deformation gradient.
             * \param &d2FdStress2: The second derivative of the flow direction w.r.t. the stress. This 
             *     is useful if one is using this expression as the flow potential and wants the jacobian of the flow direction \f$\frac{ \partial G }{\partial Stress_{IJ} }\f$
             * \param &d2FdStressdPrecedingF: The second derivative of the flow direction w.r.t. the stress and the preceding deformation gradient
             * \param tol: The tolerance used to prevent nans in the Jacobians
             */
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH(  computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatVector dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableType pressure;
            variableVector deviatoricReferenceStress;

            variableVector dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;

            variableVector d2DevStressdStressdRCG, d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            variableVector dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim );
            variableVector dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG, dRCGdPrecedingF, 1, sot_dim, sot_dim, sot_dim );

            variableVector d2DevStressdStressdPrecedingF( sot_dim * sot_dim * sot_dim, 0 );

            variableVector d2PressuredStressdPrecedingF( sot_dim * sot_dim, 0 );

            for ( unsigned int I = 0; I < deviatoricReferenceStress.size( ); I++ ){

                for ( unsigned int J = 0; J < deviatoricReferenceStress.size( ); J++ ){

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < rightCauchyGreen.size( ); L++ ){

                            d2DevStressdStressdPrecedingF[ sot_dim * sot_dim * I + sot_dim * J + K ]
                                += d2DevStressdStressdRCG[ sot_dim * sot_dim * I + sot_dim * J + L ] * dRCGdPrecedingF[ sot_dim * L + K ];

                        }

                        d2PressuredStressdPrecedingF[ sot_dim * I + J ]
                            += d2PressuredStressdRCG[ sot_dim * I + K ] * dRCGdPrecedingF[ sot_dim * K + J ];

                    }

                }

            }

            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Evaluate the jacobians
            variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

            dFdStress = tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdStress, 1, sot_dim, sot_dim, sot_dim )
                      + BAngle * dPressuredStress;

            dFdc = - AAngle;

            dFdPrecedingF = tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdPrecedingF, 1, sot_dim, sot_dim, sot_dim )
                          + BAngle * dPressuredPrecedingF;

            //Evaluate the second-order jacobians
            constantVector EYE( sot_dim * sot_dim );
            for ( unsigned int i = 0; i < sot_dim; i++ ){ EYE[ sot_dim * i + i ] = 1; }
            variableVector dDevStressDirectiondDevStress = ( EYE - tardigradeVectorTools::matrixMultiply( devStressDirection, devStressDirection, sot_dim, 1, 1, sot_dim ) ) / ( normDevStress + tol );

            d2FdStress2 = tardigradeVectorTools::matrixMultiply( dDevStressdStress, tardigradeVectorTools::matrixMultiply( dDevStressDirectiondDevStress, dDevStressdStress, sot_dim, sot_dim, sot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim, true, false );

            d2FdStressdPrecedingF = tardigradeVectorTools::matrixMultiply( tardigradeVectorTools::matrixMultiply( dDevStressDirectiondDevStress, dDevStressdStress, sot_dim, sot_dim, sot_dim, sot_dim ), dDevStressdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim, true, false )
                                  + BAngle * d2PressuredStressdPrecedingF;

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            for ( unsigned int A = 0; A < dim; A++ ){
                                for ( unsigned int B = 0; B < dim; B++ ){
                                    d2FdStressdPrecedingF[ dim * sot_dim * I + sot_dim * J + dim * K + L ] += devStressDirection[ dim * A + B ] * d2DevStressdStressdPrecedingF[ dim * sot_dim * sot_dim * A + sot_dim * sot_dim * B + dim * dim * dim * I + dim * dim * J + dim * K + L ];
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
             * \param &stressMeasure: The stress measure
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The preceding deformation gradient
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
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
                                                               variableVector &yieldValue, variableVector &dFdStress, variableVector &dFdc,
                                                               variableVector &dFdPrecedingF ){
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
             * \param &stressMeasure: The stress measure
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The preceding deformation gradient
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
             * \param &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
             * \param &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
             * \param &dFdPrecedingF: The Jacobian of the yield function w.r.t. the preceding deformation gradient
             *     tensor.
             */

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatVector dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableVector pressure;
            variableVector deviatoricReferenceStress;

            variableVector dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );

            floatVector dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, tot_dim, sot_dim, sot_dim, sot_dim );
            floatVector dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG,  dRCGdPrecedingF, dim, sot_dim, sot_dim, sot_dim );

            //Compute the l2norm of the deviatoric stress
            variableVector normDevStress;
            variableVector dNormDevStressdDevStress;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress ) );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Construct the Jacobians
            dFdStress = tardigradeVectorTools::matrixMultiply( dNormDevStressdDevStress, dDevStressdStress, dim, tot_dim, tot_dim, tot_dim )
                      + BAngle * dPressuredStress;

            dFdc = variableVector( cohesion.size( ) * cohesion.size( ), 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ dFdc[ dim * i + i ] = 1; }
            dFdc *= -AAngle;

            dFdPrecedingF = tardigradeVectorTools::matrixMultiply( dNormDevStressdDevStress, dDevStressdPrecedingF, dim, tot_dim, tot_dim, sot_dim )
                          + BAngle * dPressuredPrecedingF;

        }

        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                           const variableVector &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableVector &dFdStress, variableVector &dFdc,
                                                           variableVector &dFdPrecedingF, variableVector &d2FdStress2,
                                                           variableVector &d2FdStressdPrecedingF ){
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
             * \param &stressMeasure: The stress measure
             * \param &cohesion: The cohesion measure.
             * \param &precedingDeformationGradient: The preceding deformation gradient
             * \param &frictionAngle: The friction angle
             * \param &beta: The beta parameter
             * \param &yieldValue: The yield value.
             * \param &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
             * \param &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
             * \param &dFdPrecedingF: The Jacobian of the yield function w.r.t. the preceding deformation gradient
             * \param &d2FdStress2: The second order Jacobian of the yield function w.r.t. the reference 
             *     higher order stress.
             * \param &d2FdStressdPrecedingF: The second order Jacobian of the yield function w.r.t. the 
             *     reference higher order stress and the preceding deformation gradient
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatVector dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableVector pressure;
            variableVector deviatoricReferenceStress;

            variableVector dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;

            variableVector d2DevStressdStressdRCG, d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            floatVector dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, tot_dim, sot_dim, sot_dim, sot_dim );
            floatVector dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG,  dRCGdPrecedingF, dim, sot_dim, sot_dim, sot_dim );

            variableVector d2DevStressdStressdPrecedingF = tardigradeVectorTools::matrixMultiply( d2DevStressdStressdRCG, dRCGdPrecedingF, tot_dim * tot_dim, sot_dim, sot_dim, sot_dim );

            variableVector d2PressuredStressdPrecedingF = tardigradeVectorTools::matrixMultiply( d2PressuredStressdRCG, dRCGdPrecedingF, dim * tot_dim, sot_dim, sot_dim, sot_dim );

            //Compute the l2norm of the deviatoric stress
            variableVector normDevStress;
            variableVector dNormDevStressdDevStress;
            variableVector d2NormDevStressdDevStress2;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
                                                                                                                  dNormDevStressdDevStress,
                                                                                                                  d2NormDevStressdDevStress2 ) )

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Construct the Jacobians
            dFdStress = tardigradeVectorTools::matrixMultiply( dNormDevStressdDevStress, dDevStressdStress, dim, tot_dim, tot_dim, tot_dim )
                      + BAngle * dPressuredStress;
    
            dFdc = variableVector( cohesion.size( ) * cohesion.size( ), 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ dFdc[ dim * i + i ] = 1; }
            dFdc *= -AAngle;
    
            dFdPrecedingF = tardigradeVectorTools::matrixMultiply( dNormDevStressdDevStress, dDevStressdPrecedingF, dim, tot_dim, tot_dim, sot_dim )
                          + BAngle * dPressuredPrecedingF;
    
            //Construct the second-order jacobians
            d2FdStress2 = variableVector( dim * tot_dim * tot_dim, 0 );
            d2FdStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dNormDevStressdDevStress, d2DevStressdStressdPrecedingF, dim, tot_dim, tot_dim, tot_dim * sot_dim )
                                  + BAngle * d2PressuredStressdPrecedingF;

            for ( unsigned int K = 0; K < 3; K++ ){ //TODO: This should be able to be reduced into a matrix multiplication
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
                                                            d2FdStressdPrecedingF[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                                                += d2NormDevStressdDevStress2[ dim * dim * dim * dim * dim * dim * K + dim * dim * dim * dim * dim * Q + dim * dim * dim * dim * A + dim * dim * dim * B + dim * dim * C + dim * D + E ]
                                                                 * dDevStressdStress[ dim * dim * dim * dim * dim * Q + dim * dim * dim * dim * A + dim * dim * dim * B + dim * dim * L + dim * M + N ]
                                                                 * dDevStressdPrecedingF[ dim * dim * dim * dim * C + dim * dim * dim * D + dim * dim * E + dim * O + P ];
                                                            for ( unsigned int F = 0; F < 3; F++ ){
                                                                d2FdStress2[ dim * dim * dim * dim * dim * dim * K + dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                                                    += d2NormDevStressdDevStress2[ dim * dim * dim * dim * dim * dim * K + dim * dim * dim * dim * dim * A + dim * dim * dim * dim * B + dim * dim * dim * C + dim * dim * D + dim * E + F ]
                                                                     * dDevStressdStress[ dim * dim * dim * dim * dim * A + dim * dim * dim * dim * B + dim * dim * dim * C + dim * dim * L + dim * M + N ]
                                                                     * dDevStressdStress[ dim * dim * dim * dim * dim * D + dim * dim * dim * dim * E + dim * dim * dim * F + dim * dim * O + dim * P + Q ];
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
            constexpr unsigned int dim = 3;

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
                                                  variableVector &dPlasticMacroLdElasticRCG,
                                                  variableVector &dPlasticMacroLdMacroFlowDirection,
                                                  variableVector &dPlasticMacroLdMicroFlowDirection ){
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
             * \param &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     macro plastic multiplier.
             * \param &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     micro plastic multiplier.
             * \param &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     elastic right Cauchy-Green deformation tensor.
             * \param &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     macro flow direction tensor.
             * \param &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     micro flow direction tensor.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                     macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                     dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );
            )

            dPlasticMacroLdElasticRCG = variableVector( sot_dim * sot_dim, 0 );
            dPlasticMacroLdMacroFlowDirection = variableVector( sot_dim * sot_dim, 0 );
            dPlasticMacroLdMicroFlowDirection = variableVector( sot_dim * sot_dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                        dPlasticMacroLdMacroFlowDirection[ dim * sot_dim * Bb + sot_dim * Kb + dim * Kb + Ob ]
                            += macroGamma * inverseElasticRightCauchyGreen[ dim * Bb + Ob ];

                        dPlasticMacroLdMicroFlowDirection[ dim * sot_dim * Bb + sot_dim * Kb + dim * Kb + Ob ]
                            += microGamma * inverseElasticRightCauchyGreen[ dim * Bb + Ob ];

                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMacroLdElasticRCG[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
                                -= inverseElasticRightCauchyGreen[ dim * Bb + Ob ]
                                 * plasticMacroVelocityGradient[ dim * Pb + Kb ];

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
             * \param &microGamma: The micro plastic multiplier.
             * \param &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * \param &microFlowDirection: The micro plastic flow direction.
             * \param &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticMicroRightCauchyGreen.size() == sot_dim, "The elastic micro right Cauchy-Green deformation tensor is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticPsi.size() == sot_dim, "The elastic micro deformation tensor Psi is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( inverseElasticPsi.size() == sot_dim, "The inverse of the elastic micro deformation tensor Psi is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( microFlowDirection.size() == dim * dim, "The micro flow direction of the elastic micro plastic flow direction is not 3D" );

            plasticMicroVelocityGradient = variableVector( sot_dim, 0 );

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
             * \param &microGamma: The micro plastic multiplier.
             * \param &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * \param &microFlowDirection: The micro plastic flow direction.
             * \param &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             * \param &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro plastic multiplier.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;

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
                                                  variableVector &dPlasticMicroLdElasticMicroRCG,
                                                  variableVector &dPlasticMicroLdElasticPsi,
                                                  variableVector &dPlasticMicroLdMicroFlowDirection ){
            /*!
             * Compute the plastic micro velocity gradient
             *
             *  \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
             *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
             *
             * \param &microGamma: The micro plastic multiplier.
             * \param &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * \param &microFlowDirection: The micro plastic flow direction.
             * \param &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             * \param &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro plastic multiplier.
             * \param &dPlasticMicroLdElasticMicroRCG: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro right Cauchy-Green deformation tensor.
             * \param &dPlasticMicroLdElasticPsi: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro deformation measure Psi.
             * \param &dPlasticMicroLdMicroFlowDirection: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro flow direction.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(

                computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen,
                                                     elasticPsi, inverseElasticPsi, microFlowDirection,
                                                     plasticMicroVelocityGradient, dPlasticMicroLdMicroGamma );

            )

            //Assemble the Jacobians
            dPlasticMicroLdElasticMicroRCG = variableVector( sot_dim * sot_dim, 0 );

            dPlasticMicroLdElasticPsi = variableVector( sot_dim * sot_dim, 0 );

            dPlasticMicroLdMicroFlowDirection = variableVector( sot_dim * sot_dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMicroLdElasticPsi[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
                                -= inverseElasticPsi[ dim * Bb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Kb ];

                            for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                                dPlasticMicroLdMicroFlowDirection[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
                                    += microGamma * inverseElasticPsi[ dim * Bb + Pb ]
                                     * inverseElasticPsi[ dim * Lb + Ob ]
                                     * elasticMicroRightCauchyGreen[ dim * Lb + Kb ];

                                dPlasticMicroLdElasticMicroRCG[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Kb ]
                                    += microGamma * inverseElasticPsi[ dim * Bb + Pb ] * microFlowDirection[ dim * Lb + Pb ]
                                     * inverseElasticPsi[ dim * Ob + Lb ];

                                for ( unsigned int Eb = 0; Eb < dim; Eb++ ){

                                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                                        dPlasticMicroLdElasticPsi[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
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

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient ){
            /*!
             * Compute the plastic micro gradient velocity gradient.
             *
             * \f$\bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]\f$
             *
             * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
             * 
             * \param &microGradientGamma: The micro gradient plastic multiplier.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
             * \param &elasticGamma: The elastic higher order deformation measure Gamma.
             * \param &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
             * \param &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
             * \param &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
             */

            variableVector skewTerm;
            return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                elasticGamma, microGradientFlowDirection,
                                                                plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                                skewTerm );

        }

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm ){
            /*!
             * Compute the plastic micro gradient velocity gradient.
             *
             * \f$\bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]\f$
             *
             * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
             * 
             * \param &microGradientGamma: The micro gradient plastic multiplier.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
             * \param &elasticGamma: The elastic higher order deformation measure Gamma.
             * \param &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
             * \param &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
             * \param &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
             * \param &skewTerm: The skew term ( times 2 ) from the higher order computation.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( microGradientGamma.size() != dim ){
                    throw std::runtime_error( "The micro gradient plastic multiplier must have a length of 3" );
                }

            )

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( elasticPsi.size() != dim  * dim ){
                    throw std::runtime_error( "The elastic micro deformation measure Psi must be 3D" );
                }

            )

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( inverseElasticPsi.size() != dim  * dim ){
                    throw std::runtime_error( "The inverse elastic micro deformation measure Psi must be 3D" );
                }

            )

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( elasticGamma.size() != dim * dim * dim ){
                    throw std::runtime_error( "The elastic higher order deformation measure Gamma must be 3D" );
                }

            )

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( microGradientFlowDirection.size() != dim * dim * dim * dim ){
                    throw std::runtime_error( "The micro gradient flow direction must be 3D" );
                }

            )

            TARDIGRADE_ERROR_TOOLS_CATCH(

                if ( plasticMicroVelocityGradient.size() != dim * dim ){
                    throw std::runtime_error( "The plastic micro velocity gradient must be 3D" );
                }

            )

            //Assemble the 'skew' term
            skewTerm = variableVector( dim * dim * dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Cb = 0; Cb < dim; Cb++ ){

                            for ( unsigned int Fb = 0; Fb < dim; Fb++ ){

                                skewTerm[ dim * dim * Db + dim * Mb + Kb ]
                                    += plasticMicroVelocityGradient[ dim * Db + Cb ]
                                     * inverseElasticPsi[ dim * Cb + Fb ]
                                     * elasticGamma[ dim * dim * Fb + dim * Mb + Kb ]
                                     - plasticMicroVelocityGradient[ dim * Cb + Mb ]
                                     * inverseElasticPsi[ dim * Db + Fb ]
                                     * elasticGamma[ dim * dim * Fb + dim * Cb + Kb ];

                            }

                        }

                    }

                }

            }

            plasticMicroGradientVelocityGradient = variableVector( dim * dim * dim, 0 );

            for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                            for ( unsigned int Ib = 0; Ib < dim; Ib++ ){

                                plasticMicroGradientVelocityGradient[ dim * dim * Nb + dim * Mb + Kb ]
                                    += inverseElasticPsi[ dim * Nb + Lb ]
                                     * ( microGradientGamma[ Ib ] * microGradientFlowDirection[ dim * dim * dim * Ib + dim * dim * Kb + dim * Lb + Mb ]
                                     +   elasticPsi[ dim * Lb + Ib ] * skewTerm[ dim * dim * Ib + dim * Mb + Kb ] );

                            }

                        }

                    }

                }

            }

        }

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableVector &dPlasticMicroGradientLdPlasticMicroL ){
            /*!
             * Compute the plastic micro gradient velocity gradient.
             *
             * \f$ \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right] \f$
             *
             * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
             * 
             * \param &microGradientGamma: The micro gradient plastic multiplier.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
             * \param &elasticGamma: The elastic higher order deformation measure Gamma.
             * \param &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
             * \param &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
             * \param &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
             * \param &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the micro gradient gamma.
             * \param &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the platic micro velocity gradient.
             */

            variableVector skewTerm;
            return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi, elasticGamma,
                                                                microGradientFlowDirection, plasticMicroVelocityGradient,
                                                                plasticMicroGradientVelocityGradient, skewTerm,
                                                                dPlasticMicroGradientLdMicroGradientGamma,
                                                                dPlasticMicroGradientLdPlasticMicroL );
        }

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableVector &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableVector &dPlasticMicroGradientLdPlasticMicroL ){
            /*!
             * Compute the plastic micro gradient velocity gradient.
             *
             * \f$\bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]\f$
             *
             * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
             * 
             * \param &microGradientGamma: The micro gradient plastic multiplier.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
             * \param &elasticGamma: The elastic higher order deformation measure Gamma.
             * \param &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
             * \param &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
             * \param &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
             * \param &skewTerm: Two times the skew term.
             * \param &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the micro gradient gamma.
             * \param &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the platic micro velocity gradient.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(

                computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection,
                                                             plasticMicroVelocityGradient,
                                                             plasticMicroGradientVelocityGradient, skewTerm )

            )

            dPlasticMicroGradientLdPlasticMicroL = variableVector( tot_dim * sot_dim, 0 );
            dPlasticMicroGradientLdMicroGradientGamma = variableVector( tot_dim * dim, 0 );

            for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                            for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                                dPlasticMicroGradientLdMicroGradientGamma[ dim * dim * dim * Lb + dim * dim * Mb + dim * Kb + Ob ]
                                    += inverseElasticPsi[ dim * Lb + Pb ]
                                     * microGradientFlowDirection[ dim * dim * dim * Ob + dim * dim * Kb + dim * Pb + Mb ];

                                dPlasticMicroGradientLdPlasticMicroL[ dim * dim * sot_dim * Lb + dim * sot_dim * Mb + sot_dim * Kb + dim * Lb + Ob ]
                                    += inverseElasticPsi[ dim * Ob + Pb ] * elasticGamma[ dim * dim * Pb + dim * Mb + Kb ];

                                dPlasticMicroGradientLdPlasticMicroL[ dim * dim * sot_dim * Lb + dim * sot_dim * Mb + sot_dim * Kb + dim * Ob + Mb ]
                                    -= inverseElasticPsi[ dim * Lb + Pb ] * elasticGamma[ dim * dim * Pb + dim * Ob + Kb ];

                            }

                        }

                    }

                }

            }

        }

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableVector &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableVector &dPlasticMicroGradientLdElasticPsi,
                                                          variableVector &dPlasticMicroGradientLdElasticGamma,
                                                          variableVector &dPlasticMicroGradientLdMicroGradientFlowDirection ){
            /*!
             * Compute the plastic micro gradient velocity gradient.
             *
             * \f$\bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]\f$
             *
             * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
             * 
             * \param &microGradientGamma: The micro gradient plastic multiplier.
             * \param &elasticPsi: The elastic micro deformation measure Psi.
             * \param &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
             * \param &elasticGamma: The elastic higher order deformation measure Gamma.
             * \param &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
             * \param &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
             * \param &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
             * \param &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the micro gradient gamma.
             * \param &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
             *     velocity gradient w.r.t. the platic micro velocity gradient.
             * \param &dPlasticMicroGradientLdElasticPsi: The Jacobian of the plastic micro gradient
             *     velocity gradient w.r.t. the elastic micro deformation tensor Psi.
             * \param &dPlasticMicroGradientLdElasticGamma: The Jacobian of the plastic micro gradient
             *     velocity gradient w.r.t. the elastic higher ordrer deformation tensor Gamma.
             * \param &dPlasticMicroGradientLdMicroGradientFlowDirection: The Jacobian of the plastic micro gradient
             *     velocity gradient w.r.t. the micro gradient flow direction.
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            variableVector skewTerm;
            TARDIGRADE_ERROR_TOOLS_CATCH(

                computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection,
                                                             plasticMicroVelocityGradient,
                                                             plasticMicroGradientVelocityGradient, skewTerm,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL );

            )

            dPlasticMicroGradientLdElasticPsi = variableVector( tot_dim * sot_dim, 0 );

            dPlasticMicroGradientLdElasticGamma = variableVector( tot_dim * tot_dim, 0 );

            dPlasticMicroGradientLdMicroGradientFlowDirection = variableVector( tot_dim * dim * tot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                            for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                                dPlasticMicroGradientLdElasticPsi[ dim * dim * sot_dim * Db + dim * sot_dim * Mb + sot_dim * Kb + dim * Ob + Pb ]
                                    += -inverseElasticPsi[ dim * Db + Ob ]
                                     * plasticMicroGradientVelocityGradient[ dim * dim * Pb + dim * Mb + Kb ]
                                     + inverseElasticPsi[ dim * Db + Ob ] * skewTerm[ dim * dim * Pb + dim * Mb + Kb ];

                                dPlasticMicroGradientLdElasticGamma[ dim * dim * tot_dim * Db + dim * tot_dim * Mb + tot_dim * Kb + dim * dim * Ob + dim * Pb + Kb ]
                                    -= plasticMicroVelocityGradient[ dim * Pb + Mb ]
                                     * inverseElasticPsi[ dim * Db + Ob ];

                                dPlasticMicroGradientLdElasticGamma[ dim * dim * tot_dim * Db + dim * tot_dim * Mb + tot_dim * Kb + dim * dim * Ob + dim * Mb + Kb ]
                                    += plasticMicroVelocityGradient[ dim * Db + Pb ] * inverseElasticPsi[ dim * Pb + Ob ];

                                dPlasticMicroGradientLdMicroGradientFlowDirection[ dim * dim * tot_dim * dim * Db + dim * tot_dim * dim * Mb + tot_dim * dim * Kb + dim * dim * dim * Ob + dim * dim * Kb + dim * Pb + Mb ]
                                    += inverseElasticPsi[ dim * Db + Pb ] * microGradientGamma[ Ob ];

                                for ( unsigned int Qb = 0; Qb < dim; Qb++ ){

                                    for ( unsigned int Rb = 0; Rb < dim; Rb++ ){

                                        dPlasticMicroGradientLdElasticPsi[ dim * dim * sot_dim * Db + dim * sot_dim * Mb + sot_dim * Kb +dim * Ob + Pb ]
                                            -= plasticMicroVelocityGradient[ dim * Db + Qb ]
                                             * inverseElasticPsi[ dim * Qb + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                             * elasticGamma[ dim * dim * Rb + dim * Mb + Kb ]
                                             - plasticMicroVelocityGradient[ dim * Qb + Mb ]
                                             * inverseElasticPsi[ dim * Db + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                             * elasticGamma[ dim * dim * Rb + dim * Qb + Kb ];

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        const parameterType alpha ){
            /*!
             * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
             *
             * \param &Dt: The change in time.
             * \param &currentPlasticMicroDeformation: The inverse of the current micro deformation.
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
             *     velocity gradient.
             * \param &previousPlasticMicroDeformation: The plastic micro deformation 
             *     from the last converged increment.
             * \param &previousPlasticMicroGradient: The micro gradient deformation in the 
             *     intermediate configuation from the last converged increment.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
             *     velocity gradient from the last converged increment.
             * \param &currentPlasticMicroGradient: The current plastic micro gradient 
             *    deformation in the intermediate configuration.
             * \param alpha: The integration parameter (0 is explicit, 1 is implicit).
             */

            variableVector LHS;

            evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                       currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                       previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                       previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                       previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient, LHS,
                                       alpha );
        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableVector &LHS,
                                        const parameterType alpha ){
            /*!
             * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
             *
             * \param &Dt: The change in time.
             * \param &currentPlasticMicroDeformation: The inverse of the current micro deformation.
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
             *     velocity gradient.
             * \param &previousPlasticMicroDeformation: The the plastic micro deformation 
             *     from the last converged increment.
             * \param &previousPlasticMicroGradient: The micro gradient deformation in the 
             *     intermediate configuation from the last converged increment.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
             *     velocity gradient from the last converged increment.
             * \param &currentPlasticMicroGradient: The current plastic micro gradient 
             *    deformation in the intermediate configuration.
             * \param &LHS: The left-hand-side matrix.
             * \param alpha: The integration parameter (0 is explicit, 1 is implicit).
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
            constexpr unsigned int fot_dim = tot_dim * dim;

            TARDIGRADE_ERROR_TOOLS_CHECK( currentPlasticMicroDeformation.size() == sot_dim, "The plastic micro-deformation must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( currentPlasticMacroVelocityGradient.size() == sot_dim, "The plastic macro velocity gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( currentPlasticMicroVelocityGradient.size() == sot_dim, "The plastic micro velocity gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( currentPlasticMicroGradientVelocityGradient.size() == tot_dim, "The plastic micro gradient velocity gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( previousPlasticMicroDeformation.size() == sot_dim, "The previous plastic micro-deformation must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( previousPlasticMicroGradient.size() == tot_dim, "The previous plastic micro gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( previousPlasticMacroVelocityGradient.size() == sot_dim, "The previous plastic macro velocity gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( previousPlasticMicroVelocityGradient.size() == sot_dim, "The previous plastic micro velocity gradient must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( previousPlasticMicroGradientVelocityGradient.size() == tot_dim, "The previous plastic micro gradient velocity gradient must be 3D" );

            //Assemble the A term ( forcing term ) and the fourth order A term
            variableVector DtAtilde( tot_dim, 0 );
            variableVector previousFourthA( fot_dim, 0 );
            variableVector currentFourthA( fot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int B = 0; B < dim; B++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        previousFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Kb ]
                            += previousPlasticMicroVelocityGradient[ dim * Db + B ];

                        previousFourthA[ dim * dim * dim * Db + dim * dim * Db + dim * B + Kb ]
                            -= previousPlasticMacroVelocityGradient[ dim * Kb + B ];

                        currentFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Kb ]
                            += currentPlasticMicroVelocityGradient[ dim * Db + B ];

                        currentFourthA[ dim * dim * dim * Db + dim * dim * Db + dim * B + Kb ]
                            -= currentPlasticMacroVelocityGradient[ dim * Kb + B ];

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                            DtAtilde[ dim * dim * Db + dim * B + Kb ] += Dt
                                * ( ( 1 - alpha ) * previousPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                          *  previousPlasticMicroDeformation[ dim * Lb + B ]
                                + alpha * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                        * currentPlasticMicroDeformation[ dim * Lb + B ] );

                        }

                    }

                }

            }

            //Assemble the right-hand side and left-hand side term
            variableVector RHS = DtAtilde;
            LHS = variableVector( tot_dim * tot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                        RHS[ dim * dim * Db + dim * B + Kb ]
                           += previousPlasticMicroGradient[ dim * dim * Db + dim * B + Kb ];
                        LHS[ dim * dim * tot_dim * Db + dim * tot_dim * B + tot_dim * Kb + dim * dim * Db + dim * B + Kb ] += 1;
                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                           for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
                              RHS[ dim * dim * Db + dim * B + Kb ]
                                 += Dt * ( 1. - alpha ) * previousFourthA[ dim * dim * dim * Db + dim * dim * Bb + dim * Kb + Lb ]
                                  * previousPlasticMicroGradient[ dim * dim * Bb + dim * B + Lb ];
                              LHS[ dim * dim * tot_dim * Db + dim * tot_dim * B + tot_dim * Kb + dim * dim * Lb + dim * B + Bb ]
                                  -= Dt * alpha * currentFourthA[ dim * dim * dim * Db + dim * dim * Lb + dim * Kb + Bb ];
                           }
                        }
                    }
                }
            }

            //Solve for the current plastic micro gradient
            unsigned int rank;
            currentPlasticMicroGradient = tardigradeVectorTools::solveLinearSystem( LHS, RHS, RHS.size( ), RHS.size( ), rank );

            TARDIGRADE_ERROR_TOOLS_CHECK( rank == RHS.size(), "The left hand side matrix is not full rank" );

        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha ){
            /*!
             * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
             *
             * \param &Dt: The change in time.
             * \param &currentPlasticMicroDeformation: The inverse of the current micro deformation.
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
             *     velocity gradient.
             * \param &previousPlasticMicroDeformation: The the plastic micro deformation 
             *     from the last converged increment.
             * \param &previousPlasticMicroGradient: The micro gradient deformation in the 
             *     intermediate configuation from the last converged increment.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
             *     velocity gradient from the last converged increment.
             * \param &currentPlasticMicroGradient: The current plastic micro gradient 
             *    deformation in the intermediate configuration.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroDeformation: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro deformation.
             * \param &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic macro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro gradient velocity gradient.
             * \param alpha: The integration parameter (0 is explicit, 1 is implicit).
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
            constexpr unsigned int fot_dim = tot_dim * dim;

            //Compute the new currentPlasticMicroGradient
            variableVector LHS;
            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           LHS, alpha );
            )

            //Compute the negative partial derivatives w.r.t. currentFourthA and the current part of DtAtilde
            //We do this in vector form so that we can interface with Eigen easier
            variableVector negdRdCurrentDtAtilde( dim * dim * dim * dim * dim * dim, 0 );
            variableVector negdRdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim, 0 );

            //Also assemble jacobians of the A terms
            variableVector dCurrentDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            variableVector dCurrentDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );
            variableVector dCurrentFourthAdMacroVelocityGradient( fot_dim * sot_dim, 0 );
            variableVector dCurrentFourthAdMicroVelocityGradient( fot_dim * sot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                        negdRdCurrentDtAtilde[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Db + dim * B + Kb ] += 1;

                        dCurrentFourthAdMacroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Kb + B ] -= 1;

                        dCurrentFourthAdMicroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * B + dim * sot_dim * Kb + sot_dim * Kb + dim * Db + B ] += 1;

                        for ( unsigned int Rb = 0; Rb < dim; Rb++ ){
                            dCurrentDTAtildedPlasticMicroDeformation[ dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Rb + B ]
                                += Dt * alpha * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Rb + Kb ];

                            dCurrentDTAtildedPlasticMicroGradientVelocityGradient[ dim * dim * tot_dim * Db + dim * tot_dim * B + tot_dim * Kb + dim * dim * Db + dim * Rb + Kb ]
                                += Dt * alpha * currentPlasticMicroDeformation[ dim * Rb + B ];

                            for ( unsigned int S = 0; S < dim; S++ ){

                                negdRdCurrentFourthA[ dim * dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * dim * B + dim * dim * dim * dim * Kb + dim * dim * dim * Db + dim * dim * Rb + dim * Kb + S ]
                                    += Dt * alpha * currentPlasticMicroGradient[ dim * dim * Rb + dim * B + S ];
                            }
                        }
                    }
                }
            }

            //Solve for the Jacobians
            variableVector dCurrentPlasticMicroGradientdCurrentDTAtilde( dim * dim * dim * dim * dim * dim );
            variableVector dCurrentPlasticMicroGradientdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim );

            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > LHSMat( LHS.data(), tot_dim, tot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCDA( negdRdCurrentDtAtilde.data(), tot_dim, tot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCFA( negdRdCurrentFourthA.data(), tot_dim, sot_dim * sot_dim );

            Eigen::ColPivHouseholderQR< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > qrSolver( LHSMat );

            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X1( dCurrentPlasticMicroGradientdCurrentDTAtilde.data(), tot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X2( dCurrentPlasticMicroGradientdCurrentFourthA.data(), tot_dim, sot_dim * sot_dim );

            X1 = qrSolver.solve( nDRDCDA );
            X2 = qrSolver.solve( nDRDCFA );

            //Assemble the final terms of the deformation
            dCurrentPlasticMicroGradientdPlasticMicroDeformation = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                          dCurrentDTAtildedPlasticMicroDeformation, tot_dim, tot_dim, tot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                                               dCurrentFourthAdMacroVelocityGradient, tot_dim, sot_dim * sot_dim, sot_dim * sot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                                               dCurrentFourthAdMicroVelocityGradient, tot_dim, sot_dim * sot_dim, sot_dim * sot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                                       dCurrentDTAtildedPlasticMicroGradientVelocityGradient, tot_dim, tot_dim, tot_dim, tot_dim );

        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                        variableVector &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                                        variableVector &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha ){
            /*!
             * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
             *
             * \param &Dt: The change in time.
             * \param &currentPlasticMicroDeformation: The inverse of the current micro deformation.
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
             *     velocity gradient.
             * \param &previousPlasticMicroDeformation: The the plastic micro deformation 
             *     from the last converged increment.
             * \param &previousPlasticMicroGradient: The micro gradient deformation in the 
             *     intermediate configuation from the last converged increment.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
             *     from the last converged increment.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
             *     velocity gradient from the last converged increment.
             * \param &currentPlasticMicroGradient: The current plastic micro gradient 
             *    deformation in the intermediate configuration.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroDeformation: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro deformation.
             * \param &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic macro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the plastic micro gradient velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation: The jacobian of the plastic 
             *     micro deformation w.r.t. the previous plastic micro deformation.
             * \param &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the intermediate configuration spatial gradient of the plastic micro deformation from the last
             *     converged increment
             * \param &dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the previous plastic macro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the previous plastic micro velocity gradient.
             * \param &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient: The jacobian of the plastic 
             *     micro deformation w.r.t. the previous plastic micro gradient velocity gradient.
             * \param alpha: The integration parameter (0 is explicit, 1 is implicit).
             */

            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
            constexpr unsigned int fot_dim = tot_dim * dim;

            //Compute the required identity terms
            constantVector eye( sot_dim, 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ eye[ dim * i + i ] = 1; }

            //Compute the new currentPlasticMicroGradient
            variableVector LHS;
            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           LHS, alpha );
            )

            //Compute the negative partial derivatives w.r.t. currentFourthA and the current part of DtAtilde
            //We do this in vector form so that we can interface with Eigen easier
            variableVector negdRdCurrentDtAtilde( tot_dim * tot_dim, 0 );
            variableVector negdRdCurrentFourthA( tot_dim * fot_dim, 0 );

            //Also assemble jacobians of the A terms
            variableVector dCurrentDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            variableVector dCurrentDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );
            variableVector dCurrentFourthAdMacroVelocityGradient( fot_dim * sot_dim, 0 );
            variableVector dCurrentFourthAdMicroVelocityGradient( fot_dim * sot_dim, 0 );

            variableVector dPreviousDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            variableVector dPreviousDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );

            variableVector dRHSdPreviousPlasticMicroGradient( tot_dim * tot_dim, 0 );
            variableVector dRHSdPreviousPlasticMacroVelocityGradient( tot_dim * sot_dim, 0 );
            variableVector dRHSdPreviousPlasticMicroVelocityGradient( tot_dim * sot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                        negdRdCurrentDtAtilde[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Db + dim * B + Kb ] += 1;

                        dCurrentFourthAdMicroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * B + dim * sot_dim * Kb + sot_dim * Kb + dim * Db + B ] += 1;

                        dCurrentFourthAdMacroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Kb + B ] -= 1;

                        for ( unsigned int Rb = 0; Rb < dim; Rb++ ){

                            dCurrentDTAtildedPlasticMicroDeformation[ dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Rb + B ]
                                += Dt * alpha * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Rb + Kb ];

                            dPreviousDTAtildedPlasticMicroDeformation[ dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Rb + B ]
                                += Dt * ( 1. - alpha ) * previousPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Rb + Kb ];

                            dRHSdPreviousPlasticMacroVelocityGradient[ dim * dim * dim * dim * Db + dim * dim * dim * B + dim * dim * Kb + dim * Rb + Kb ]
                                -= Dt * ( 1. - alpha ) * previousPlasticMicroGradient[ dim * dim * Db + dim * B + Rb ];

                            dRHSdPreviousPlasticMicroVelocityGradient[ dim * dim * dim * dim * Db + dim * dim * dim * B + dim * dim * Kb + dim * Db + Rb ]
                                += Dt * ( 1. - alpha ) * previousPlasticMicroGradient[ dim * dim * Rb + dim * B + Kb ];

                            dCurrentDTAtildedPlasticMicroGradientVelocityGradient[ dim * dim * tot_dim * Db + dim * tot_dim * B + tot_dim * Kb + dim * dim * Db + dim * Rb + Kb ]
                                += Dt * alpha * currentPlasticMicroDeformation[ dim * Rb + B ];

                            dPreviousDTAtildedPlasticMicroGradientVelocityGradient[ dim * dim * tot_dim * Db + dim * tot_dim * B + tot_dim * Kb + dim * dim * Db + dim * Rb + Kb ]
                                += Dt * ( 1. - alpha ) * previousPlasticMicroDeformation[ dim * Rb + B ];

                            for ( unsigned int S = 0; S < dim; S++ ){

                                negdRdCurrentFourthA[ dim * dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * dim * B + dim * dim * dim * dim * Kb + dim * dim * dim * Db + dim * dim * Rb + dim * Kb + S ]
                                    += Dt * alpha * currentPlasticMicroGradient[ dim * dim * Rb + dim * B + S ];

                                for ( unsigned int Tb = 0; Tb < dim; Tb++ ){
                                    dRHSdPreviousPlasticMicroGradient[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Rb + dim * S + Tb ]
                                        += ( eye[ dim * Db + Rb ] * eye[ dim * Kb + Tb ] * eye[ dim * B + S ] + Dt * ( 1 - alpha ) * ( previousPlasticMicroVelocityGradient[ dim * Db + Rb ] * eye[ dim * Kb + Tb ]
                                                                                                              -   previousPlasticMacroVelocityGradient[ dim * Tb + Kb ] * eye[ dim * Db + Rb ] ) * eye[ dim * B + S ] );

                                }
                            }
                        }
                    }
                }
            }

            floatVector dCurrentPlasticMicroGradientdCurrentDTAtilde( tot_dim * tot_dim, 0 );
            floatVector dCurrentPlasticMicroGradientdCurrentFourthA( tot_dim * fot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient = floatVector( tot_dim * tot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient = floatVector( tot_dim * sot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient = floatVector( tot_dim * sot_dim, 0 );

            //Solve for the Jacobians
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > LHSMat( LHS.data(), tot_dim, tot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCDA( negdRdCurrentDtAtilde.data(), tot_dim, tot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCFA( negdRdCurrentFourthA.data(), tot_dim, fot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DRDPPMG( dRHSdPreviousPlasticMicroGradient.data(), tot_dim, tot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DRDPPMaVG( dRHSdPreviousPlasticMacroVelocityGradient.data(), tot_dim, sot_dim );
            Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DRDPPMiVG( dRHSdPreviousPlasticMicroVelocityGradient.data(), tot_dim, sot_dim );

            Eigen::ColPivHouseholderQR< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > qrSolver( LHSMat );

            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X1( dCurrentPlasticMicroGradientdCurrentDTAtilde.data(), tot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X2( dCurrentPlasticMicroGradientdCurrentFourthA.data(), tot_dim, fot_dim );
            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DCPMGDPMG( dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient.data(), tot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DCPMGDPMaVG( dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient.data(), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > DCPMGDPMiVG( dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient.data(), tot_dim, sot_dim );

            X1          = qrSolver.solve( nDRDCDA );
            X2          = qrSolver.solve( nDRDCFA );
            DCPMGDPMG   = qrSolver.solve( DRDPPMG );
            DCPMGDPMaVG = qrSolver.solve( DRDPPMaVG );
            DCPMGDPMiVG = qrSolver.solve( DRDPPMiVG );

            //Assemble the final terms of the deformation
            dCurrentPlasticMicroGradientdPlasticMicroDeformation = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                          dCurrentDTAtildedPlasticMicroDeformation, tot_dim, tot_dim, tot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                                               dCurrentFourthAdMacroVelocityGradient, tot_dim, fot_dim, fot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                                               dCurrentFourthAdMicroVelocityGradient, tot_dim, fot_dim, fot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                                       dCurrentDTAtildedPlasticMicroGradientVelocityGradient,
                                                                                                                       tot_dim, tot_dim, tot_dim, tot_dim );

            dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                                  dPreviousDTAtildedPlasticMicroDeformation,
                                                                                                                  tot_dim, tot_dim, tot_dim, sot_dim );

            dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient = tardigradeVectorTools::matrixMultiply( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                                                               dPreviousDTAtildedPlasticMicroGradientVelocityGradient,
                                                                                                                               tot_dim, tot_dim, tot_dim, tot_dim );

        }

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
            /*!
             * Evolve the plastic deformation
             *
             * \param &Dt: The timestep
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
             *     velocity gradient.
             * \param &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
             *     converged timestep.
             * \param &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
             *     timestep.
             * \param &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
             *     timestep.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
             *     last converged timestep.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
             *     last converged timestep.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
             *     at the end of the last converged timestep.
             * \param &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
             * \param &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
             * \param &currentPlasticMicroGradient: The current value of the plastic micro gradient.
             * \param alphaMacro: The integration parameter for the macro plasticity. 0 explicit, 1 implicit. Defaults to 0.5.
             * \param alphaMicro: The integration parameter for the micro plasticity. 0 explicit, 1 implicit. Defaults to 0.5.
             * \param alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
             */

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                            currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                            1. - alphaMacro, 1 );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            1. - alphaMicro, 1 );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           alphaMicroGradient );
            )

        }

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableVector &dPlasticFdPlasticMacroL,
                                       variableVector &dPlasticMicroDeformationdPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPlasticMacroL,
                                       variableVector &dPlasticMicroGradientdPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPlasticMicroGradientL,
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
            /*!
             * Evolve the plastic deformation
             *
             * \param &Dt: The timestep
             * \param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * \param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * \param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
             *     velocity gradient.
             * \param &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
             *     converged timestep.
             * \param &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
             *     timestep.
             * \param &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
             *     timestep.
             * \param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
             *     last converged timestep.
             * \param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
             *     last converged timestep.
             * \param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
             *     at the end of the last converged timestep.
             * \param &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
             * \param &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
             * \param &currentPlasticMicroGradient: The current value of the plastic micro gradient.
             * \param &dPlasticFdPlasticMacroL: The Jacobian of the plastic deformation gradient w.r.t. the plastic 
             *     macro velocity gradient.
             * \param &dPlasticMicroDeformationdPlasticMicroL: The Jacobian of the plastic micro-deformation w.r.t. 
             *     the plastic micro velocity gradient.
             * \param &dPlasticMicroGradientdPlasticMacroL: The Jacobian of the plastic micro gradient deformation 
             *     w.r.t. the plastic macro velocity gradient.
             * \param &dPlasticMicroGradientdPlasticMicroL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the plastic micro velocity gradient.
             * \param &dPlasticMicroGradientdPlasticMicroGradientL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the plastic micro gradient velocity gradient.
             * \param alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
             * \param alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
             * \param alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
             */

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                            currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                            dPlasticFdPlasticMacroL, 1. - alphaMacro, 1 );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, 1. - alphaMicro, 1 );
            )

            variableVector dPlasticMicroGradientdPlasticMicroDeformation;
            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           dPlasticMicroGradientdPlasticMicroDeformation,
                                           dPlasticMicroGradientdPlasticMacroL, dPlasticMicroGradientdPlasticMicroL,
                                           dPlasticMicroGradientdPlasticMicroGradientL, alphaMicroGradient );
            )

            dPlasticMicroGradientdPlasticMicroL += tardigradeVectorTools::matrixMultiply( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                                          dPlasticMicroDeformationdPlasticMicroL, tot_dim, sot_dim, sot_dim, sot_dim );

        }

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableVector &dPlasticFdPlasticMacroL,
                                       variableVector &dPlasticMicroDeformationdPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPlasticMacroL,
                                       variableVector &dPlasticMicroGradientdPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPlasticMicroGradientL,
                                       variableVector &dPlasticFdPreviousPlasticF,
                                       variableVector &dPlasticFdPreviousPlasticMacroL,
                                       variableVector &dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                       variableVector &dPlasticMicroDeformationdPreviousPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                       variableVector &dPlasticMicroGradientdPreviousPlasticMicroGradient,
                                       variableVector &dPlasticMicroGradientdPreviousPlasticMacroL,
                                       variableVector &dPlasticMicroGradientdPreviousPlasticMicroL,
                                       variableVector &dPlasticMicroGradientdPreviousPlasticMicroGradientL,
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
            /*!
             * Evolve the plastic deformation
             *
             * :param &Dt: The timestep
             * :param &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
             * :param &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
             * :param &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
             *     velocity gradient.
             * :param &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
             *     converged timestep.
             * :param &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
             *     timestep.
             * :param &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
             *     timestep.
             * :param &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
             *     last converged timestep.
             * :param &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
             *     last converged timestep.
             * :param &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
             *     at the end of the last converged timestep.
             * :param &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
             * :param &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
             * :param &currentPlasticMicroGradient: The current value of the plastic micro gradient.
             * :param &dPlasticFdPlasticMacroL: The Jacobian of the plastic deformation gradient w.r.t. the plastic 
             *     macro velocity gradient.
             * :param &dPlasticMicroDeformationdPlasticMicroL: The Jacobian of the plastic micro-deformation w.r.t. 
             *     the plastic micro velocity gradient.
             * :param &dPlasticMicroGradientdPlasticMacroL: The Jacobian of the plastic micro gradient deformation 
             *     w.r.t. the plastic macro velocity gradient.
             * :param &dPlasticMicroGradientdPlasticMicroL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the plastic micro velocity gradient.
             * :param &dPlasticMicroGradientdPlasticMicroGradientL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the plastic micro gradient velocity gradient.
             * :param &dPlasticFdPreviousPlasticF: The Jacobian of the plastic deformation gradient w.r.t. the previous
             *     plastic deformation gradient.
             * :param &dPlasticFdPreviousPlasticMacroL: The Jacobian of the plastic deformation gradient w.r.t. the previous plastic 
             *     macro velocity gradient.
             * :param &dPlasticMicroDeformationdPreviousPlasticMicroDeformation: The Jacobian of the plastic micro deformation w.r.t. the previous
             *     plastic micro deformation
             * :param &dPlasticMicroDeformationdPreviousPlasticMicroL: The Jacobian of the plastic micro-deformation w.r.t. 
             *     the previous plastic micro velocity gradient.
             * :param &dPlasticMicroGradientdPreviousPlasticMicroDeformation: The Jacobian of the plastic micro gradient deformation 
             *     w.r.t. the previous plastic micro deformation.
             * :param &dPlasticMicroGradientdPreviousPlasticMicroGradient: The Jacobian of the plastic micro gradient deformation 
             *     w.r.t. the previous spatial gradient in the intermediate configuration of the plastic macro deformation.
             * :param &dPlasticMicroGradientdPreviousPlasticMacroL: The Jacobian of the plastic micro gradient deformation 
             *     w.r.t. the previous plastic macro velocity gradient.
             * :param &dPlasticMicroGradientdPreviousPlasticMicroL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the previous plastic micro velocity gradient.
             * :param &dPlasticMicroGradientdPreviousPlasticMicroGradientL: The Jacobian of the plastic micro gradient deformation
             *     w.r.t. the previous plastic micro gradient velocity gradient.
             * :param alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
             * :param alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
             * :param alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
             */

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                            currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                            dPlasticFdPlasticMacroL, dPlasticFdPreviousPlasticF, dPlasticFdPreviousPlasticMacroL, 1. - alphaMacro, 1 );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPreviousPlasticMicroL, 1. - alphaMicro, 1 );
            )

            variableVector dPlasticMicroGradientdPlasticMicroDeformation;
            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           dPlasticMicroGradientdPlasticMicroDeformation,
                                           dPlasticMicroGradientdPlasticMacroL, dPlasticMicroGradientdPlasticMicroL,
                                           dPlasticMicroGradientdPlasticMicroGradientL,
                                           dPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                           dPlasticMicroGradientdPreviousPlasticMicroGradient,
                                           dPlasticMicroGradientdPreviousPlasticMacroL,
                                           dPlasticMicroGradientdPreviousPlasticMicroL,
                                           dPlasticMicroGradientdPreviousPlasticMicroGradientL,
                                           alphaMicroGradient );
            )

            dPlasticMicroGradientdPlasticMicroL += tardigradeVectorTools::matrixMultiply( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                                          dPlasticMicroDeformationdPlasticMicroL, tot_dim, sot_dim, sot_dim, sot_dim );

            dPlasticMicroGradientdPreviousPlasticMicroDeformation += tardigradeVectorTools::matrixMultiply( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                                                            dPlasticMicroDeformationdPreviousPlasticMicroDeformation, tot_dim, sot_dim, sot_dim, sot_dim );

            dPlasticMicroGradientdPreviousPlasticMicroL += tardigradeVectorTools::matrixMultiply( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                                                  dPlasticMicroDeformationdPreviousPlasticMicroL, tot_dim, sot_dim, sot_dim, sot_dim );


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

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

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
            floatVector PK2Stress(                     stress->begin( ),               stress->begin( ) + 1 * sot_dim );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * sot_dim, stress->begin( ) + 2 * sot_dim );

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * sot_dim, stress->begin( ) + 2 * sot_dim + tot_dim );

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

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *stress;

            floatVector dFpdSubFs;

            const floatVector *dF1dF;

            const floatVector *dF1dFn;

            floatVector dChipdSubChis;

            const floatVector *dChi1dChi;

            const floatVector *dChi1dChin;

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
            floatVector dFpdF(  sot_dim * sot_dim, 0 );

            floatVector dFpdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            floatVector dChipdChi(  sot_dim * sot_dim, 0 );

            floatVector dChipdChin( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dFpdF[ sot_dim * i + j ] += dFpdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dF )[ sot_dim * k + j ];

                        dChipdChi[ sot_dim * i + j ] += dChipdSubChis[ num_configs * sot_dim * i +  k ] * ( *dChi1dChi )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: Investigate ways to prevent cache misses here

                    dFpdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dFpdSubFs[ num_configs * sot_dim * i + j + sot_dim ];

                    dChipdChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dChipdSubChis[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dFpdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dFpdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dFn )[ ( num_configs - 1 ) * sot_dim * k + j ];

                        dChipdChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dChipdSubChis[ num_configs * sot_dim * i + k ] * ( *dChi1dChin )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

            }

            // Extract the stresses from the stress vector
            floatVector PK2Stress(                     stress->begin( ),               stress->begin( ) + 1 * sot_dim );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * sot_dim, stress->begin( ) + 2 * sot_dim );;

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * sot_dim, stress->begin( ) + 2 * sot_dim + tot_dim );

            // Push the stresses forward to the current configuration of the plastic configuration
            floatVector macroDrivingStress;

            floatVector symmetricMicroDrivingStress;

            floatVector higherOrderDrivingStress;

            floatVector dMacrodFp;

            floatVector dMacrodPK2;

            floatVector dMicrodFp;

            floatVector dMicrodSigma;

            floatVector dHigherdFp;

            floatVector dHigherdChip;

            floatVector dHigherdM;

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

                set_previousdMacroDrivingStressdF( tardigradeVectorTools::matrixMultiply( dMacrodFp, dFpdF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdMacroDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dMacrodFp, dFpdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_previousdSymmetricMicroDrivingStressdF( tardigradeVectorTools::matrixMultiply( dMicrodFp, dFpdF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdSymmetricMicroDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dMicrodFp, dFpdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_previousdHigherOrderDrivingStressdF( tardigradeVectorTools::matrixMultiply( dHigherdFp, dFpdF, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdHigherOrderDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dHigherdFp, dFpdFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdHigherOrderDrivingStressdChi( tardigradeVectorTools::matrixMultiply( dHigherdChip, dChipdChi, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdHigherOrderDrivingStressdChin( tardigradeVectorTools::matrixMultiply( dHigherdChip, dChipdChin, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

                set_dMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_dMacroDrivingStressdF( tardigradeVectorTools::matrixMultiply( dMacrodFp, dFpdF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dMacroDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dMacrodFp, dFpdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_dSymmetricMicroDrivingStressdF( tardigradeVectorTools::matrixMultiply( dMicrodFp, dFpdF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dSymmetricMicroDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dMicrodFp, dFpdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_dHigherOrderDrivingStressdF( tardigradeVectorTools::matrixMultiply( dHigherdFp, dFpdF, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dHigherOrderDrivingStressdFn( tardigradeVectorTools::matrixMultiply( dHigherdFp, dFpdFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dHigherOrderDrivingStressdChi( tardigradeVectorTools::matrixMultiply( dHigherdChip, dChipdChi, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dHigherOrderDrivingStressdChin( tardigradeVectorTools::matrixMultiply( dHigherdChip, dChipdChin, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

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

            floatVector dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

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

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            const floatVector *dMacroDrivingStressdStress;

            const floatVector *dMacroDrivingStressdF;

            const floatVector *dMacroDrivingStressdFn;

            const floatVector *dMicroDrivingStressdStress;

            const floatVector *dMicroDrivingStressdF;

            const floatVector *dMicroDrivingStressdFn;

            const floatVector *dMicroGradientDrivingStressdStress;

            const floatVector *dMicroGradientDrivingStressdF;

            const floatVector *dMicroGradientDrivingStressdFn;

            const floatVector *dMicroGradientDrivingStressdChi;

            const floatVector *dMicroGradientDrivingStressdChin;

            const floatVector *precedingDeformationGradient;

            const floatVector *dPrecedingFdF;

            const floatVector *dPrecedingFdFn;

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

            floatVector dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

            floatVector d2MacroFlowdDrivingStress2,         d2MacroFlowdDrivingStressdPrecedingF,
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

            floatVector d2MacroFlowdDrivingStressdMacroStress( sot_dim * sot_dim, 0 );

            floatVector d2MicroFlowdDrivingStressdMicroStress( sot_dim * sot_dim, 0 );

            floatVector d2MicroGradientFlowdDrivingStressdMicroGradientStress( dim * tot_dim * tot_dim, 0 );

            floatVector d2MacroFlowdDrivingStressdF( sot_dim * sot_dim, 0 );

            floatVector d2MicroFlowdDrivingStressdF( sot_dim * sot_dim, 0 );

            floatVector d2MicroGradientFlowdDrivingStressdF( dim * tot_dim * sot_dim, 0 );

            floatVector d2MicroGradientFlowdDrivingStressdChi( dim * tot_dim * sot_dim, 0 );

            floatVector d2MacroFlowdDrivingStressdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            floatVector d2MicroFlowdDrivingStressdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            floatVector d2MicroGradientFlowdDrivingStressdFn( dim * tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            floatVector d2MicroGradientFlowdDrivingStressdChin( dim * tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int I = 0; I < sot_dim; I++ ){

                for ( unsigned int K = 0; K < sot_dim; K++ ){

                    for ( unsigned int J = 0; J < sot_dim; J++ ){

                        d2MacroFlowdDrivingStressdMacroStress[ sot_dim * I + J ]
                            += d2MacroFlowdDrivingStress2[ sot_dim * I + K ] * ( *dMacroDrivingStressdStress )[ sot_dim * K + J ];

                        d2MicroFlowdDrivingStressdMicroStress[ sot_dim * I  + J ]
                            += d2MicroFlowdDrivingStress2[ sot_dim * I + K ] * ( *dMicroDrivingStressdStress )[ sot_dim * K + J ];

                        d2MacroFlowdDrivingStressdF[ sot_dim * I + J ]
                            += d2MacroFlowdDrivingStress2[ sot_dim * I + K ] * ( *dMacroDrivingStressdF )[ sot_dim * K + J ]
                             + d2MacroFlowdDrivingStressdPrecedingF[ sot_dim * I + K ] * ( *dPrecedingFdF )[ sot_dim * K + J ];

                        d2MicroFlowdDrivingStressdF[ sot_dim * I + J ]
                            += d2MicroFlowdDrivingStress2[ sot_dim * I + K ] * ( *dMicroDrivingStressdF )[ sot_dim * K + J ]
                             + d2MicroFlowdDrivingStressdPrecedingF[ sot_dim * I + K ] * ( *dPrecedingFdF )[ sot_dim * K + J ];

                    }

                }

                for ( unsigned int J = 0; J < sot_dim; J++ ){

                    for ( unsigned int K = 0; K < ( num_configs - 1 ) * sot_dim; K++ ){

                        d2MacroFlowdDrivingStressdFn[ ( num_configs - 1 ) * sot_dim * I + K ]
                            += d2MacroFlowdDrivingStressdPrecedingF[ sot_dim * I + J ] * ( *dPrecedingFdFn )[ ( num_configs - 1 ) * sot_dim * J + K ]
                             + d2MacroFlowdDrivingStress2[ sot_dim * I + J ] * ( *dMacroDrivingStressdFn )[ ( num_configs - 1 ) * sot_dim * J + K ];

                        d2MicroFlowdDrivingStressdFn[ ( num_configs - 1 ) * sot_dim * I + K ]
                            += d2MicroFlowdDrivingStressdPrecedingF[ sot_dim * I + J ] * ( *dPrecedingFdFn )[ ( num_configs - 1 ) * sot_dim * J + K ]
                             + d2MicroFlowdDrivingStress2[ sot_dim * I + J ] * ( *dMicroDrivingStressdFn ) [ ( num_configs - 1 ) * sot_dim * J + K ];

                    }

                }

            }

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int J = 0; J < tot_dim; J++ ){

                    for ( unsigned int L = 0; L < tot_dim; L++ ){

                        for ( unsigned int K = 0; K < tot_dim; K++ ){

                            d2MicroGradientFlowdDrivingStressdMicroGradientStress[ tot_dim * tot_dim * I + tot_dim * J + K ]
                                += d2MicroGradientFlowdDrivingStress2[ tot_dim * tot_dim * I + tot_dim * J + L ] * ( *dMicroGradientDrivingStressdStress )[ tot_dim * L + K ];

                        }

                    }
                    for ( unsigned int K = 0; K < tot_dim; K++ ){

                        for ( unsigned int L = 0; L < sot_dim; L++ ){

                            d2MicroGradientFlowdDrivingStressdF[ tot_dim * sot_dim * I + sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ tot_dim * tot_dim * I + tot_dim * J + K ] * ( *dMicroGradientDrivingStressdF )[ sot_dim * K + L ];

                            d2MicroGradientFlowdDrivingStressdChi[ tot_dim * sot_dim * I + sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ tot_dim * tot_dim * I + tot_dim * J + K ] * ( *dMicroGradientDrivingStressdChi )[ sot_dim * K + L ];

                        }

                    }
                    for ( unsigned int K = 0; K < tot_dim; K++ ){

                        for ( unsigned int L = 0; L < ( num_configs - 1 ) * sot_dim; L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ tot_dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ tot_dim * tot_dim * I + tot_dim * J + K ] * ( *dMicroGradientDrivingStressdFn ) [ ( num_configs - 1 ) * sot_dim * K + L ];

                            d2MicroGradientFlowdDrivingStressdChin[ tot_dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ tot_dim * tot_dim * I + tot_dim * J + K ] * ( *dMicroGradientDrivingStressdChin ) [ ( num_configs - 1 ) * sot_dim * K + L ];

                        }

                    }

                    for ( unsigned int K = 0; K < sot_dim; K++ ){

                        for ( unsigned int L = 0; L < sot_dim; L++ ){

                            d2MicroGradientFlowdDrivingStressdF[ tot_dim * sot_dim * I + sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStressdPrecedingF[ tot_dim * sot_dim * I + sot_dim * J + K ] * ( *dPrecedingFdF )[ sot_dim * K + L ];

                        }

                    }

                    for ( unsigned int K = 0; K < sot_dim; K++ ){

                        for ( unsigned int L = 0; L < ( num_configs - 1 ) * sot_dim; L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ tot_dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + L ]
                                += d2MicroGradientFlowdDrivingStressdPrecedingF[ tot_dim * sot_dim * I + sot_dim * J + K ] * ( *dPrecedingFdFn )[ ( num_configs - 1 ) * sot_dim * K + L ];

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

            const unsigned int num_pisvs = get_plasticStateVariables( )->size( );

            const floatVector *plasticStrainLikeISVs;

            if ( isPrevious ){

                plasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            }
            else{

                plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            }

            const unsigned int num_psisvs = plasticStrainLikeISVs->size( );

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

            floatVector dMacroCohesiondISVs( num_pisvs, 0 );

            floatVector dMicroCohesiondISVs( num_pisvs, 0 );

            floatVector dMicroGradientCohesiondISVs( ( get_plasticStrainLikeISVs( )->size( ) - 2 ) * num_pisvs, 0 );

            floatType macroCohesion           = ( *get_macroHardeningParameters( ) )[ 0 ] + ( *get_macroHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 0 ];

            dMacroCohesiondISVs[ get_plasticMultipliers( )->size( ) + 0 ] = ( *get_macroHardeningParameters( ) )[ 1 ];

            floatType microCohesion           = ( *get_microHardeningParameters( ) )[ 0 ] + ( *get_microHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 1 ];

            dMicroCohesiondISVs[ get_plasticMultipliers( )->size( ) + 1 ] = ( *get_microHardeningParameters( ) )[ 1 ];

            floatVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * floatVector( plasticStrainLikeISVs->begin( ) + 2,
                                                                                                                                                                     plasticStrainLikeISVs->end( ) );

            for ( unsigned int i = 2; i < num_psisvs; i++ ){

                dMicroGradientCohesiondISVs[ num_pisvs * ( i - 2 ) + get_plasticMultipliers( )->size( ) + i ] = ( *get_microGradientHardeningParameters( ) )[ 1 ];

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

            const floatVector *dMicroGradientFlowdc;

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

            const unsigned int num_pm = plasticMultipliers->size( );

            floatVector evolutionRates( plasticMultipliers->size( ), 0 );

            evolutionRates[ 0 ] = -( *dMacroFlowdc ) * ( *plasticMultipliers )[ 0 ];

            evolutionRates[ 1 ] = -( *dMicroFlowdc ) * ( *plasticMultipliers )[ 1 ];

            for ( unsigned int i = 2; i < num_pm; i++ ){

                for ( unsigned int j = 2; j < num_pm; j++ ){

                    evolutionRates[ i ] -= ( *dMicroGradientFlowdc )[ 3 * ( i - 2 ) + j - 2 ] * ( *plasticMultipliers )[ j ];

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

            const floatVector *dMicroGradientFlowdc;

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

            const unsigned int num_pm = plasticMultipliers->size( );

            const unsigned int num_psvs = get_plasticStateVariables( )->size( );

            floatVector evolutionRates( num_pm, 0 );

            floatVector dEvolutionRatesdStateVariables( num_pm * num_psvs, 0 );

            evolutionRates[ 0 ] = -( *dMacroFlowdc ) * ( *plasticMultipliers )[ 0 ];

            evolutionRates[ 1 ] = -( *dMicroFlowdc ) * ( *plasticMultipliers )[ 1 ];

            dEvolutionRatesdStateVariables[ num_psvs * 0 + 0 ] = -( *dMacroFlowdc );

            dEvolutionRatesdStateVariables[ num_psvs * 1 + 1 ] = -( *dMicroFlowdc );

            for ( unsigned int i = 2; i < num_pm; i++ ){

                for ( unsigned int j = 2; j < num_pm; j++ ){

                    evolutionRates[ i ] -= ( *dMicroGradientFlowdc )[ 3 * ( i - 2 ) + j - 2 ] * ( *plasticMultipliers )[ j ];

                    dEvolutionRatesdStateVariables[ num_psvs * i + j ] = -( *dMicroGradientFlowdc )[ 3 * ( i - 2 ) + j - 2 ];

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

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, 1 - ( *getIntegrationParameter( ) ) ) );

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

            const unsigned int num_isvs = previousPlasticStrainLikeISVs->size( );
            const unsigned int num_psvs = get_plasticStateVariables( )->size( );

            floatVector dISVs, updatedISVs;

            floatVector dISVsdEvolutionRates;

            if ( addPrevious ){

                floatVector dISVsdPreviousEvolutionRates;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolutionFlatJ( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, dISVsdEvolutionRates, dISVsdPreviousEvolutionRates, 1 - ( *getIntegrationParameter( ) ) ) );

                floatVector dISVsdStateVariables = tardigradeVectorTools::matrixMultiply( dISVsdPreviousEvolutionRates, *get_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( ),
                                                                                          num_isvs, num_isvs, num_isvs, num_psvs );

                for ( unsigned int i = 0; i < updatedISVs.size( ); i++ ){

                    dISVsdStateVariables[ num_psvs * i  + i + ( *getNumPlasticMultipliers( ) ) ] += 1;

                }

                set_dUpdatedPlasticStrainLikeISVsdPreviousStateVariables( dISVsdStateVariables );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolutionFlatJ( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, dISVsdEvolutionRates, 1 - ( *getIntegrationParameter( ) ) ) );

            }

            set_dUpdatedPlasticStrainLikeISVsdStateVariables( tardigradeVectorTools::matrixMultiply( dISVsdEvolutionRates, *get_dPlasticStrainLikeISVEvolutionRatesdStateVariables( ),
                                                                                                     num_isvs, num_isvs, num_isvs, num_psvs ) );

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

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *dMacroCohesiondStateVariables;

            const floatVector *dMacroDrivingStressdStress;

            const floatVector *dMacroDrivingStressdF;

            const floatVector *dMacroDrivingStressdFn;

            const floatVector *dMicroCohesiondStateVariables;

            const floatVector *dMicroDrivingStressdStress;

            const floatVector *dMicroDrivingStressdF;

            const floatVector *dMicroDrivingStressdFn;

            const floatVector *dMicroGradientCohesiondStateVariables;

            const floatVector *dMicroGradientDrivingStressdStress;

            const floatVector *dMicroGradientDrivingStressdF;

            const floatVector *dMicroGradientDrivingStressdFn;

            const floatVector *dMicroGradientDrivingStressdChi;

            const floatVector *dMicroGradientDrivingStressdChin;

            const floatVector *precedingDeformationGradient;

            const floatVector *dPrecedingFdF;

            const floatVector *dPrecedingFdFn;

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

            floatVector dMicroGradientYielddDrivingStress;

            floatVector dMicroGradientYielddCohesion;

            floatVector dMicroGradientYielddPrecedingF;

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

            floatVector dMacroYielddStress = tardigradeVectorTools::matrixMultiply( dMacroYielddDrivingStress, *dMacroDrivingStressdStress, 1, sot_dim, sot_dim, sot_dim );

            floatVector dMacroYielddStateVariables = dMacroYielddCohesion * ( *dMacroCohesiondStateVariables );

            floatVector dMacroYielddF = tardigradeVectorTools::matrixMultiply( dMacroYielddDrivingStress, *dMacroDrivingStressdF, 1, sot_dim, sot_dim, sot_dim )
                                      + tardigradeVectorTools::matrixMultiply( dMacroYielddPrecedingF,    *dPrecedingFdF,         1, sot_dim, sot_dim, sot_dim);

            floatVector dMacroYielddFn = tardigradeVectorTools::matrixMultiply( dMacroYielddDrivingStress, *dMacroDrivingStressdFn, 1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                       + tardigradeVectorTools::matrixMultiply( dMacroYielddPrecedingF,    *dPrecedingFdFn,         1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dMicroYielddStress = tardigradeVectorTools::matrixMultiply( dMicroYielddDrivingStress, *dMicroDrivingStressdStress, 1, sot_dim, sot_dim, sot_dim );

            floatVector dMicroYielddStateVariables = dMicroYielddCohesion * ( *dMicroCohesiondStateVariables );

            floatVector dMicroYielddF = tardigradeVectorTools::matrixMultiply( dMicroYielddDrivingStress, *dMicroDrivingStressdF, 1, sot_dim, sot_dim, sot_dim )
                                      + tardigradeVectorTools::matrixMultiply( dMicroYielddPrecedingF,    *dPrecedingFdF,         1, sot_dim, sot_dim, sot_dim );

            floatVector dMicroYielddFn = tardigradeVectorTools::matrixMultiply( dMicroYielddDrivingStress, *dMicroDrivingStressdFn, 1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                       + tardigradeVectorTools::matrixMultiply( dMicroYielddPrecedingF,    *dPrecedingFdFn,         1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dMicroGradientYielddStress = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdStress, 3, tot_dim, tot_dim, tot_dim );

            floatVector dMicroGradientYielddStateVariables = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddCohesion, *dMicroGradientCohesiondStateVariables, 3, 3, 3, dMicroGradientCohesiondStateVariables->size( ) / 3 );

            floatVector dMicroGradientYielddF = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdF, 3, tot_dim, tot_dim, sot_dim )
                                              + tardigradeVectorTools::matrixMultiply( dMicroGradientYielddPrecedingF, *dPrecedingFdF, 3, sot_dim, sot_dim, sot_dim );

            floatVector dMicroGradientYielddFn = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdFn, 3, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                                               + tardigradeVectorTools::matrixMultiply( dMicroGradientYielddPrecedingF, *dPrecedingFdFn, 3, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dMicroGradientYielddChi = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChi, 3, tot_dim, tot_dim, sot_dim );

            floatVector dMicroGradientYielddChin = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChin, 3, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

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

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector dPrecedingFdSubFs;

            const floatVector *dF1dF;

            const floatVector *dF1dFn;

            if ( isPrevious ){

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF  = hydra->get_previousdF1dF( );

                dF1dFn = hydra->get_previousdF1dFn( );

            }
            else{

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF  = hydra->get_dF1dF( );

                dF1dFn = hydra->get_dF1dFn( );

            }

            // Construct the derivatives of the preceding F

            floatVector dPrecedingFdF( sot_dim * sot_dim, 0 );

            floatVector dPrecedingFdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dPrecedingFdF[ sot_dim * i + j ] += dPrecedingFdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dF )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: See if this can be sped up by accessing the cache better

                    dPrecedingFdFn[ ( num_configs - 1 ) * sot_dim * i + j ] = dPrecedingFdSubFs[ num_configs * sot_dim * i + sot_dim + j ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dPrecedingFdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dPrecedingFdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dFn )[ ( num_configs - 1 ) * sot_dim * k + j ];

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

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector dPrecedingChidSubChis;

            const floatVector *dChi1dChi;

            const floatVector *dChi1dChin;

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

            floatVector dPrecedingChidChi( sot_dim * sot_dim, 0 );

            floatVector dPrecedingChidChin( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dPrecedingChidChi[ sot_dim * i + j ] += dPrecedingChidSubChis[ num_configs * sot_dim * i + k ] * ( *dChi1dChi )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: See if this can be sped up by accessing the cache better

                    dPrecedingChidChin[ ( num_configs - 1 ) * sot_dim * i + j ] = dPrecedingChidSubChis[ num_configs * sot_dim * i + sot_dim + j ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dPrecedingChidChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dPrecedingChidSubChis[ num_configs * sot_dim * i + k ] * ( *dChi1dChin )[ ( num_configs - 1 ) * sot_dim * k + j ];

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

        void residual::setPrecedingGradientMicroDeformation( ){
            /*!
             * Set the preceding gradient of the micro-deformation
             */

            setPrecedingGradientMicroDeformation( false );

        }

        void residual::setPreviousPrecedingGradientMicroDeformation( ){
            /*!
             * Set the previous preceding gradient of the micro-deformation
             */

            setPrecedingGradientMicroDeformation( true );

        }

        void residual::setPrecedingGradientMicroDeformation( const bool isPrevious ){
            /*!
             * Set the preceding gradient w.r.t. the micro deformation
             * 
             * \param isPrevious: Flag for whether to set the current (flase) or previous (true) gradient of the micro deformation
             */

            const unsigned int tot_dim = hydra->getTOTDimension( );

            if ( isPrevious ){

                set_previousPrecedingGradientMicroDeformation( floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                                            hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression

            }
            else{

                set_precedingGradientMicroDeformation( floatVector( hydra->get_gradientMicroConfigurations( )->begin( ),
                                                                    hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression

            }

        }

        void residual::setdPrecedingGradientMicroDeformationdFn( ){
            /*!
             * Set the Jacobian of the preceding gradient of the micro-deformation w.r.t. the sub configurations
             */

            setPrecedingGradientMicroDeformationJacobians( false );

        }

        void residual::setdPrecedingGradientMicroDeformationdChi( ){
            /*!
             * Set the Jacobian of the preceding gradient of the micro-deformation w.r.t. the micro deformation
             */

            setPrecedingGradientMicroDeformationJacobians( false );

        }

        void residual::setdPrecedingGradientMicroDeformationdChin( ){
            /*!
             * Set the Jacobian of the preceding gradient of the micro-deformation w.r.t. the sub-micro deformations
             */

            setPrecedingGradientMicroDeformationJacobians( false );

        }

        void residual::setdPrecedingGradientMicroDeformationdGradChi( ){
            /*!
             * Set the Jacobian of the preceding gradient of the micro-deformation w.r.t. the spatial gradient of the micro deformation
             */

            setPrecedingGradientMicroDeformationJacobians( false );

        }

        void residual::setdPrecedingGradientMicroDeformationdGradChin( ){
            /*!
             * Set the Jacobian of the preceding gradient of the micro-deformation w.r.t. the spatial local gradient of the sub-micro deformations
             */

            setPrecedingGradientMicroDeformationJacobians( false );

        }

        void residual::setPreviousdPrecedingGradientMicroDeformationdFn( ){
            /*!
             * Set the previous Jacobian of the preceding gradient of the micro-deformation w.r.t. the sub configurations
             */

            setPrecedingGradientMicroDeformationJacobians( true );

        }

        void residual::setPreviousdPrecedingGradientMicroDeformationdChi( ){
            /*!
             * Set the previous Jacobian of the preceding gradient of the micro-deformation w.r.t. the micro deformation
             */

            setPrecedingGradientMicroDeformationJacobians( true );

        }

        void residual::setPreviousdPrecedingGradientMicroDeformationdChin( ){
            /*!
             * Set the previous Jacobian of the preceding gradient of the micro-deformation w.r.t. the sub-micro deformations
             */

            setPrecedingGradientMicroDeformationJacobians( true );

        }

        void residual::setPreviousdPrecedingGradientMicroDeformationdGradChi( ){
            /*!
             * Set the previous Jacobian of the preceding gradient of the micro-deformation w.r.t. the spatial gradient of the micro deformation
             */

            setPrecedingGradientMicroDeformationJacobians( true );

        }

        void residual::setPreviousdPrecedingGradientMicroDeformationdGradChin( ){
            /*!
             * Set the previous Jacobian of the preceding gradient of the micro-deformation w.r.t. the spatial local gradient of the sub-micro deformations
             */

            setPrecedingGradientMicroDeformationJacobians( true );

        }

        void residual::setPrecedingGradientMicroDeformationJacobians( const bool isPrevious ){
            /*!
             * Set the preceding gradient w.r.t. the micro deformation
             * 
             * \param isPrevious: Flag for whether to set the current (flase) or previous (true) gradient of the micro deformation
             */

            const unsigned int tot_dim = hydra->getTOTDimension( );

            if ( isPrevious ){

                set_previousPrecedingGradientMicroDeformation( floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                                            hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression

                set_previousdPrecedingGradientMicroDeformationdFn(       *hydra->get_previousdGradChi1dFn( ) );

                set_previousdPrecedingGradientMicroDeformationdChi(      *hydra->get_previousdGradChi1dChi( ) );

                set_previousdPrecedingGradientMicroDeformationdChin(     *hydra->get_previousdGradChi1dChin( ) );

                set_previousdPrecedingGradientMicroDeformationdGradChi(  *hydra->get_previousdGradChi1dGradChi( ) );

                set_previousdPrecedingGradientMicroDeformationdGradChin( *hydra->get_previousdGradChi1dGradChin( ) );

            }
            else{

                set_precedingGradientMicroDeformation( floatVector( hydra->get_gradientMicroConfigurations( )->begin( ),
                                                                    hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression
                                                                                                                                       //
                set_dPrecedingGradientMicroDeformationdFn(       *hydra->get_dGradChi1dFn( ) );

                set_dPrecedingGradientMicroDeformationdChi(      *hydra->get_dGradChi1dChi( ) );

                set_dPrecedingGradientMicroDeformationdChin(     *hydra->get_dGradChi1dChin( ) );

                set_dPrecedingGradientMicroDeformationdGradChi(  *hydra->get_dGradChi1dGradChi( ) );

                set_dPrecedingGradientMicroDeformationdGradChin( *hydra->get_dGradChi1dGradChin( ) );

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

        void residual::setPlasticGradientMicroVelocityGradient( ){
            /*!
             * Set the plastic gradient micro velocity gradient
             */

            setPlasticVelocityGradients( false );

        }

        void residual::setPreviousPlasticGradientMicroVelocityGradient( ){
            /*!
             * Set the previous plastic gradient macro velocity gradient
             */

            setPlasticVelocityGradients( true );

        }

        void residual::setPlasticVelocityGradients( const bool isPrevious ){
            /*!
             * Set the plastic macro velocity gradient
             * 
             * \param isPrevious: Flag for if the previous (true) or current (false) velocity gradient should be calculated
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const floatVector *precedingDeformationGradient;

            const floatVector *precedingMicroDeformation;

            const floatVector *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const floatVector *dMacroFlowdDrivingStress;

            const floatVector *dMicroFlowdDrivingStress;

            const floatVector *dMicroGradientFlowdDrivingStress;

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                precedingMicroDeformation = get_previousPrecedingMicroDeformation( );

                precedingGradientMicroDeformation = get_previousPrecedingGradientMicroDeformation( );

                plasticMultipliers = get_previousPlasticMultipliers( );

                dMacroFlowdDrivingStress = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_previousdMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_previousdMicroGradientFlowdDrivingStress( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

                precedingMicroDeformation = get_precedingMicroDeformation( );

                precedingGradientMicroDeformation = get_precedingGradientMicroDeformation( );

                plasticMultipliers = get_plasticMultipliers( );

                dMacroFlowdDrivingStress = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_dMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_dMicroGradientFlowdDrivingStress( );

            }

            // Form the preceding RCG and its inverse
            floatVector precedingRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG ) );

            floatVector inversePrecedingRCG = precedingRCG;
            Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( inversePrecedingRCG.data(), 3, 3 );
            mat = mat.inverse( ).eval( );

            // Form the precedingPsi and its inverse
            floatVector precedingPsi;

            TARDIGRADE_ERROR_TOOLS_CATCH( precedingPsi = tardigradeVectorTools::matrixMultiply( *precedingDeformationGradient, *precedingMicroDeformation, dim, dim, dim, dim, true, false ) );

            floatVector inversePrecedingPsi = precedingPsi;
            new (&mat ) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( inversePrecedingPsi.data(), 3, 3 );
            mat = mat.inverse( ).eval( );

            // Form the preceding micro RCG and its inverse
            floatVector precedingMicroRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG ) );

            // Form Gamma
            floatVector precedingGamma( tot_dim, 0 );
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        for ( unsigned int K = 0; K < dim; K++ ){

                            precedingGamma[ sot_dim * I + dim * J + K ]
                                += ( *precedingDeformationGradient )[ dim * i + I ] * ( *precedingGradientMicroDeformation )[ sot_dim * i + dim * J + K ];

                        }

                    }

                }

            }

            floatVector macroVelocityGradient;

            floatVector microVelocityGradient;

            floatVector gradientMicroVelocityGradient;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( ( *plasticMultipliers )[ 0 ], ( *plasticMultipliers )[ 1 ],
                                                     inversePrecedingRCG, *dMacroFlowdDrivingStress, *dMicroFlowdDrivingStress,
                                                     macroVelocityGradient )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( ( *plasticMultipliers )[ 1 ], precedingMicroRCG, precedingPsi, inversePrecedingPsi,
                                                     *dMicroFlowdDrivingStress, microVelocityGradient );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroGradientVelocityGradient( floatVector( plasticMultipliers->begin( ) + 2, plasticMultipliers->end( ) ),
                                                             precedingPsi, inversePrecedingPsi, precedingGamma,
                                                             *dMicroGradientFlowdDrivingStress,
                                                             microVelocityGradient, gradientMicroVelocityGradient );
            )

            if ( isPrevious ){

                set_previousPlasticMacroVelocityGradient( macroVelocityGradient );

                set_previousPlasticMicroVelocityGradient( microVelocityGradient );

                set_previousPlasticGradientMicroVelocityGradient( gradientMicroVelocityGradient );

            }
            else{

                set_plasticMacroVelocityGradient( macroVelocityGradient );

                set_plasticMicroVelocityGradient( microVelocityGradient );

                set_plasticGradientMicroVelocityGradient( gradientMicroVelocityGradient );

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

        void residual::setdPlasticMicroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the deformation gradient
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the deformation gradient
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

        void residual::setdPlasticMicroVelocityGradientdChi( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the micro deformation
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdChi( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the micro deformation
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

        void residual::setdPlasticMicroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the state variables
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the state variables
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdHigherOrderStress( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the higher order stress
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdHigherOrderStress( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the higher order stress
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the deformation gradient
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the deformation gradient
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdChi( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the micro deformation
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdChi( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the micro deformation
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdChin( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the sub micro deformation
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdChin( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the sub micro deformation
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdGradChi( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the spatial gradient of the micro deformation
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdGradChi( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the spatial gradient of the sub micro deformation
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdGradChin( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the local spatial gradient of the sub micro deformations
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdGradChin( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the local spatial gradient of the sub micro deformations
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticGradientMicroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the plastic gradient micro velocity gradient w.r.t. the local spatial gradient of the state variables
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticGradientMicroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the previous plastic gradient micro velocity gradient w.r.t. the local spatial gradient of the state variables
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setPlasticVelocityGradientsJacobians( const bool isPrevious ){
            /*!
             * Set the plastic macro velocity gradient and their Jacobians
             * 
             * \param isPrevious: Flag for if the previous (true) or current (false) velocity gradient should be calculated
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *precedingDeformationGradient;

            const floatVector *dPrecedingFdF;

            const floatVector *dPrecedingFdFn;

            const floatVector *dPrecedingChidChi;

            const floatVector *dPrecedingChidChin;

            const floatVector *dPrecedingGradChidFn;

            const floatVector *dPrecedingGradChidChi;

            const floatVector *dPrecedingGradChidChin;

            const floatVector *dPrecedingGradChidGradChi;

            const floatVector *dPrecedingGradChidGradChin;

            const floatVector *precedingMicroDeformation;

            const floatVector *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const floatVector *dMacroFlowdDrivingStress;

            const floatVector *dMicroFlowdDrivingStress;

            const floatVector *dMicroGradientFlowdDrivingStress;

            const floatVector *d2MacroFlowdDrivingStressdStress;

            const floatVector *d2MacroFlowdDrivingStressdF;

            const floatVector *d2MacroFlowdDrivingStressdFn;

            const floatVector *d2MicroFlowdDrivingStressdF;

            const floatVector *d2MicroFlowdDrivingStressdFn;

            const floatVector *d2MicroFlowdDrivingStressdStress;

            const floatVector *d2MicroGradientFlowdDrivingStressdStress;

            const floatVector *d2MicroGradientFlowdDrivingStressdF;

            const floatVector *d2MicroGradientFlowdDrivingStressdFn;

            const floatVector *d2MicroGradientFlowdDrivingStressdChi;

            const floatVector *d2MicroGradientFlowdDrivingStressdChin;

            if ( isPrevious ){

                dPrecedingFdF = get_previousdPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_previousdPrecedingDeformationGradientdFn( );

                dPrecedingChidChi = get_previousdPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_previousdPrecedingMicroDeformationdChin( );

                dPrecedingGradChidFn       = get_previousdPrecedingGradientMicroDeformationdFn( );

                dPrecedingGradChidChi      = get_previousdPrecedingGradientMicroDeformationdChi( );

                dPrecedingGradChidChin     = get_previousdPrecedingGradientMicroDeformationdChin( );

                dPrecedingGradChidGradChi  = get_previousdPrecedingGradientMicroDeformationdGradChi( );

                dPrecedingGradChidGradChin = get_previousdPrecedingGradientMicroDeformationdGradChin( );

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                precedingMicroDeformation = get_previousPrecedingMicroDeformation( );

                precedingGradientMicroDeformation = get_previousPrecedingGradientMicroDeformation( );

                plasticMultipliers = get_previousPlasticMultipliers( );

                d2MacroFlowdDrivingStressdStress = get_previousd2MacroFlowdDrivingStressdStress( );

                d2MacroFlowdDrivingStressdF      = get_previousd2MacroFlowdDrivingStressdF( );

                d2MacroFlowdDrivingStressdFn     = get_previousd2MacroFlowdDrivingStressdFn( );

                d2MicroFlowdDrivingStressdStress = get_previousd2MicroFlowdDrivingStressdStress( );

                d2MicroFlowdDrivingStressdF      = get_previousd2MicroFlowdDrivingStressdF( );

                d2MicroFlowdDrivingStressdFn     = get_previousd2MicroFlowdDrivingStressdFn( );

                dMacroFlowdDrivingStress         = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_previousdMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_previousdMicroGradientFlowdDrivingStress( );

                d2MicroGradientFlowdDrivingStressdStress = get_previousd2MicroGradientFlowdDrivingStressdStress( );

                d2MicroGradientFlowdDrivingStressdF      = get_previousd2MicroGradientFlowdDrivingStressdF( );

                d2MicroGradientFlowdDrivingStressdFn     = get_previousd2MicroGradientFlowdDrivingStressdFn( );

                d2MicroGradientFlowdDrivingStressdChi    = get_previousd2MicroGradientFlowdDrivingStressdChi( );

                d2MicroGradientFlowdDrivingStressdChin   = get_previousd2MicroGradientFlowdDrivingStressdChin( );

            }
            else{

                dPrecedingFdF = get_dPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_dPrecedingDeformationGradientdFn( );

                dPrecedingChidChi = get_dPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_dPrecedingMicroDeformationdChin( );

                dPrecedingGradChidFn       = get_dPrecedingGradientMicroDeformationdFn( );

                dPrecedingGradChidChi      = get_dPrecedingGradientMicroDeformationdChi( );

                dPrecedingGradChidChin     = get_dPrecedingGradientMicroDeformationdChin( );

                dPrecedingGradChidGradChi  = get_dPrecedingGradientMicroDeformationdGradChi( );

                dPrecedingGradChidGradChin = get_dPrecedingGradientMicroDeformationdGradChin( );

                precedingDeformationGradient = get_precedingDeformationGradient( );

                precedingMicroDeformation = get_precedingMicroDeformation( );

                precedingGradientMicroDeformation = get_precedingGradientMicroDeformation( );

                plasticMultipliers = get_plasticMultipliers( );

                d2MacroFlowdDrivingStressdStress = get_d2MacroFlowdDrivingStressdStress( );

                d2MacroFlowdDrivingStressdF      = get_d2MacroFlowdDrivingStressdF( );

                d2MacroFlowdDrivingStressdFn     = get_d2MacroFlowdDrivingStressdFn( );

                d2MicroFlowdDrivingStressdStress = get_d2MicroFlowdDrivingStressdStress( );

                d2MicroFlowdDrivingStressdF      = get_d2MicroFlowdDrivingStressdF( );

                d2MicroFlowdDrivingStressdFn     = get_d2MicroFlowdDrivingStressdFn( );

                dMacroFlowdDrivingStress         = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_dMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_dMicroGradientFlowdDrivingStress( );

                d2MicroGradientFlowdDrivingStressdStress = get_d2MicroGradientFlowdDrivingStressdStress( );

                d2MicroGradientFlowdDrivingStressdF      = get_d2MicroGradientFlowdDrivingStressdF( );

                d2MicroGradientFlowdDrivingStressdFn     = get_d2MicroGradientFlowdDrivingStressdFn( );

                d2MicroGradientFlowdDrivingStressdChi    = get_d2MicroGradientFlowdDrivingStressdChi( );

                d2MicroGradientFlowdDrivingStressdChin   = get_d2MicroGradientFlowdDrivingStressdChin( );

            }

            // Form the preceding RCG and its inverse
            floatVector precedingRCG;

            floatVector dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG, dRCGdPrecedingF ) );

            floatVector inversePrecedingRCG = precedingRCG;
            Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( inversePrecedingRCG.data(), 3, 3 );
            mat = mat.inverse( ).eval( );

            floatVector dRCGdF = tardigradeVectorTools::matrixMultiply( dRCGdPrecedingF,  *dPrecedingFdF, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dRCGdFn = tardigradeVectorTools::matrixMultiply( dRCGdPrecedingF, *dPrecedingFdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            // Form the precedingPsi and its inverse
            floatVector precedingPsi;

            TARDIGRADE_ERROR_TOOLS_CATCH( precedingPsi = tardigradeVectorTools::matrixMultiply( *precedingDeformationGradient, *precedingMicroDeformation, dim, dim, dim, dim, true, false ) );

            floatVector dPsidPrecedingF( sot_dim * sot_dim, 0 );

            floatVector dPsidPrecedingChi( sot_dim * sot_dim, 0 );

            floatVector eye( sot_dim, 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ eye[ dim * i + i ] = 1; }

            for ( unsigned int I = 0; I < dim; I++ ){ //TODO: Try to improve the cache access here

                for ( unsigned int J = 0; J < dim; J++ ){

                    for ( unsigned int A = 0; A < dim; A++ ){

                        for ( unsigned int B = 0; B < dim; B++ ){

                            dPsidPrecedingF[   dim * dim * dim * I + dim * dim * J + dim * A + B ] += eye[ dim * I + B ] * ( *precedingMicroDeformation )[ dim * A + J ];

                            dPsidPrecedingChi[ dim * dim * dim * I + dim * dim * J + dim * A + B ] += ( *precedingDeformationGradient )[ dim * A + I ] * eye[ dim * J + B ];

                        }

                    }

                }

            }

            floatVector dPsidF    = tardigradeVectorTools::matrixMultiply( dPsidPrecedingF,   *dPrecedingFdF     , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPsidFn   = tardigradeVectorTools::matrixMultiply( dPsidPrecedingF,   *dPrecedingFdFn    , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPsidChi  = tardigradeVectorTools::matrixMultiply( dPsidPrecedingChi, *dPrecedingChidChi , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPsidChin = tardigradeVectorTools::matrixMultiply( dPsidPrecedingChi, *dPrecedingChidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector inversePrecedingPsi = precedingPsi;
            new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( inversePrecedingPsi.data(), 3, 3 );
            mat = mat.inverse( ).eval( );

            // Form the preceding micro RCG and its inverse
            floatVector precedingMicroRCG;

            floatVector dMicroRCGdPrecedingChi;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG, dMicroRCGdPrecedingChi ) );

            floatVector dMicroRCGdChi  = tardigradeVectorTools::matrixMultiply( dMicroRCGdPrecedingChi, *dPrecedingChidChi, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dMicroRCGdChin = tardigradeVectorTools::matrixMultiply( dMicroRCGdPrecedingChi, *dPrecedingChidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            // Form Gamma
            floatVector precedingGamma( tot_dim, 0 );

            floatVector dPrecedingGammadPrecedingF( tot_dim * sot_dim, 0 );

            floatVector dPrecedingGammadPrecedingGradChi( tot_dim * tot_dim, 0 );

            for ( unsigned int i = 0; i < dim; i++ ){ //TODO: Try to improve the cache access here

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        for ( unsigned int K = 0; K < dim; K++ ){

                            precedingGamma[ sot_dim * I + dim * J + K ]
                                += ( *precedingDeformationGradient )[ dim * i + I ] * ( *precedingGradientMicroDeformation )[ sot_dim * i + dim * J + K ];

                            for ( unsigned int A = 0; A < dim; A++ ){

                                dPrecedingGammadPrecedingF[ sot_dim * sot_dim * I + dim * sot_dim * J + sot_dim * K + dim * i + A ]
                                    += eye[ dim * I + A ] * ( *precedingGradientMicroDeformation )[ sot_dim * i + dim * J + K ];

                                for ( unsigned int B = 0; B < dim; B++ ){

                                    dPrecedingGammadPrecedingGradChi[ sot_dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + sot_dim * i + dim * A + B ]
                                        += ( *precedingDeformationGradient )[ dim * i + I ] * eye[ dim * J + A ] * eye[ dim * K + B ];

                                }

                            }

                        }

                    }

                }

            }

            floatVector dPrecedingGammadF        = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingF, *dPrecedingFdF, tot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPrecedingGammadFn       = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingF,       *dPrecedingFdFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                 + tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingGradChi, *dPrecedingGradChidFn, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPrecedingGammadChi      = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingGradChi, *dPrecedingGradChidChi, tot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dPrecedingGammadChin     = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingGradChi, *dPrecedingGradChidChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPrecedingGammadGradChi  = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingGradChi, *dPrecedingGradChidGradChi, tot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dPrecedingGammadGradChin = tardigradeVectorTools::matrixMultiply( dPrecedingGammadPrecedingGradChi, *dPrecedingGradChidGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector macroVelocityGradient;

            floatVector microVelocityGradient;

            floatVector gradientMicroVelocityGradient;

            floatVector dPlasticMacroLdMacroGamma;

            floatVector dPlasticMacroLdMicroGamma;

            floatVector dPlasticMacroLdPrecedingRCG;

            floatVector dPlasticMacroLdMacroFlowDirection;

            floatVector dPlasticMacroLdMicroFlowDirection;

            floatVector dPlasticMicroLdMicroGamma;

            floatVector dPlasticMicroLdPrecedingMicroRCG;

            floatVector dPlasticMicroLdPrecedingPsi;

            floatVector dPlasticMicroLdMicroFlowDirection;

            floatVector dPlasticGradientMicroLdMicroGradientGamma;

            floatVector dPlasticGradientMicroLdPlasticMicroL;

            floatVector dPlasticGradientMicroLdPrecedingPsi;

            floatVector dPlasticGradientMicroLdPrecedingGamma;

            floatVector dPlasticGradientMicroLdMicroGradientFlowDirection;

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

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroGradientVelocityGradient( floatVector( plasticMultipliers->begin( ) + 2, plasticMultipliers->end( ) ),
                                                             precedingPsi, inversePrecedingPsi, precedingGamma,
                                                             *dMicroGradientFlowdDrivingStress,
                                                             microVelocityGradient, gradientMicroVelocityGradient,
                                                             dPlasticGradientMicroLdMicroGradientGamma,
                                                             dPlasticGradientMicroLdPlasticMicroL,
                                                             dPlasticGradientMicroLdPrecedingPsi,
                                                             dPlasticGradientMicroLdPrecedingGamma,
                                                             dPlasticGradientMicroLdMicroGradientFlowDirection );
            )

            // Assemble the Jacobians
            floatVector dPlasticMacroLdMacroStress = tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdStress, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticMacroLdMicroStress = tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdStress, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticMicroLdMicroStress = tardigradeVectorTools::matrixMultiply( dPlasticMicroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdStress, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticGradientMicroLdMicroStress = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdMicroStress, tot_dim, sot_dim, sot_dim, sot_dim );

//            floatVector reshaped_d2MicroGradientFlowdDrivingStressdStress( dim * tot_dim * tot_dim, 0 );
//
//            floatVector reshaped_d2MicroGradientFlowdDrivingStressdF( dim * tot_dim * sot_dim, 0 );
//
//            floatVector reshaped_d2MicroGradientFlowdDrivingStressdFn( dim * tot_dim * ( num_configs - 1 ) * sot_dim, 0 );
//
//            floatVector reshaped_d2MicroGradientFlowdDrivingStressdChi( dim * tot_dim * sot_dim, 0 );
//
//            floatVector reshaped_d2MicroGradientFlowdDrivingStressdChin( dim * tot_dim * ( num_configs - 1 ) * sot_dim, 0 );
//
//            for ( unsigned int i = 0; i < dim; i++ ){ //TODO: Try and improve the cache access here
//
//                for ( unsigned int j = 0; j < tot_dim; j++ ){
//
//                    for ( unsigned int k = 0; k < tot_dim; k++ ){
//
//                        reshaped_d2MicroGradientFlowdDrivingStressdStress[ tot_dim * tot_dim * i + tot_dim * j + k ]
//                            = ( *d2MicroGradientFlowdDrivingStressdStress )[ tot_dim * tot_dim * i + tot_dim * j + k ];
//
//                    }
//
//                }
//                for ( unsigned int j = 0; j < tot_dim; j++ ){
//
//                    for ( unsigned int k = 0; k < sot_dim; k++ ){
//
//                        reshaped_d2MicroGradientFlowdDrivingStressdF[ tot_dim * sot_dim * i + sot_dim * j + k ]
//                            = ( *d2MicroGradientFlowdDrivingStressdF )[ i ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) *  ( *dim ) * j + k ];
//
//                        reshaped_d2MicroGradientFlowdDrivingStressdChi[ ( *dim ) * ( *dim ) * ( *dim ) * i + j ][ k ]
//                            = ( *d2MicroGradientFlowdDrivingStressdChi )[ i ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) *  ( *dim ) * j + k ];
//
//                    }
//
//                }
//                for ( unsigned int j = 0; j < tot_dim; j++ ){
//
//                    for ( unsigned int k = 0; k < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); k++ ){
//
//                        reshaped_d2MicroGradientFlowdDrivingStressdFn[ ( *dim ) * ( *dim ) * ( *dim ) * i + j ][ k ]
//                            = ( *d2MicroGradientFlowdDrivingStressdFn )[ i ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) *  ( *dim ) * j + k ];
//
//                        reshaped_d2MicroGradientFlowdDrivingStressdChin[ ( *dim ) * ( *dim ) * ( *dim ) * i + j ][ k ]
//                            = ( *d2MicroGradientFlowdDrivingStressdChin )[ i ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) *  ( *dim ) * j + k ];
//
//                    }
//
//                }
//
//            }

            floatVector dPlasticGradientMicroLdHigherOrderStress = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdMicroGradientFlowDirection, *d2MicroGradientFlowdDrivingStressdStress,
                                                                                                          tot_dim, 3 * tot_dim, 3 * tot_dim, tot_dim );

            floatVector dPlasticMacroLdF             = tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdF,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdF,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMacroLdPrecedingRCG, dRCGdF,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticMacroLdFn            = tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdFn,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdFn,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMacroLdPrecedingRCG, dRCGdFn,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPlasticMicroLdF             = tardigradeVectorTools::matrixMultiply( dPlasticMicroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdF,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingPsi, dPsidF,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticMicroLdFn            = tardigradeVectorTools::matrixMultiply( dPlasticMicroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdFn,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingPsi, dPsidFn,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPlasticGradientMicroLdFn    = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdFn,
                                                                                              tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingPsi, dPsidFn,
                                                                                              tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadFn,
                                                                                              tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdMicroGradientFlowDirection, *d2MicroGradientFlowdDrivingStressdFn,
                                                                                              tot_dim, 3 * tot_dim, 3 * tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPlasticGradientMicroLdF     = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdF,
                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingPsi, dPsidF,
                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadF,
                                                                                              tot_dim, tot_dim, tot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdMicroGradientFlowDirection, *d2MicroGradientFlowdDrivingStressdF,
                                                                                              tot_dim, 3 * tot_dim, 3 * tot_dim, sot_dim );

            floatVector dPlasticMicroLdChi           = tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingPsi, dPsidChi,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingMicroRCG, dMicroRCGdChi,
                                                                                              sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPlasticMicroLdChin          = tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingPsi, dPsidChin,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticMicroLdPrecedingMicroRCG, dMicroRCGdChin,
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPlasticGradientMicroLdChi   = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdChi,
                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingPsi, dPsidChi,
                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadChi,
                                                                                              tot_dim, tot_dim, tot_dim, sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdMicroGradientFlowDirection, *d2MicroGradientFlowdDrivingStressdChi,
                                                                                              tot_dim, 3 * tot_dim, 3 * tot_dim, sot_dim );

            floatVector dPlasticGradientMicroLdChin  = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdChin,
                                                                                              tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingPsi, dPsidChin,
                                                                                              tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadChin,
                                                                                              tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdMicroGradientFlowDirection, *d2MicroGradientFlowdDrivingStressdChin,
                                                                                              tot_dim, 3 * tot_dim, 3 * tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPlasticGradientMicroLdGradChi  = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadGradChi,
                                                                                                 tot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dPlasticGradientMicroLdGradChin = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPrecedingGamma, dPrecedingGammadGradChin,
                                                                                                 tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            floatVector dPlasticMacroLdISVs( sot_dim * num_isvs, 0 );

            floatVector dPlasticMicroLdISVs( sot_dim * num_isvs, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                dPlasticMacroLdISVs[ num_isvs * i + 0 ] = dPlasticMacroLdMacroGamma[ i ];

                dPlasticMacroLdISVs[ num_isvs * i + 1 ] = dPlasticMacroLdMicroGamma[ i ];

                dPlasticMicroLdISVs[ num_isvs * i + 1 ] = dPlasticMicroLdMicroGamma[ i ];

            }

            floatVector dPlasticGradientMicroLdISVs = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroLdPlasticMicroL, dPlasticMicroLdISVs,
                                                                                             tot_dim, sot_dim, sot_dim, num_isvs );

            for ( unsigned int i = 0; i < tot_dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    dPlasticGradientMicroLdISVs[ num_isvs * i + j + 2 ] = dPlasticGradientMicroLdMicroGradientGamma[ 3 * i + j ];

                }

            }

            if ( isPrevious ){

                set_previousPlasticMacroVelocityGradient( macroVelocityGradient );

                set_previousPlasticMicroVelocityGradient( microVelocityGradient );

                set_previousPlasticGradientMicroVelocityGradient( gradientMicroVelocityGradient );

                set_previousdPlasticMacroVelocityGradientdMacroStress( dPlasticMacroLdMacroStress );

                set_previousdPlasticMacroVelocityGradientdMicroStress( dPlasticMacroLdMicroStress );

                set_previousdPlasticMicroVelocityGradientdMicroStress( dPlasticMicroLdMicroStress );

                set_previousdPlasticMacroVelocityGradientdF( dPlasticMacroLdF );

                set_previousdPlasticMacroVelocityGradientdFn( dPlasticMacroLdFn );

                set_previousdPlasticMicroVelocityGradientdF( dPlasticMicroLdF );

                set_previousdPlasticMicroVelocityGradientdFn( dPlasticMicroLdFn );

                set_previousdPlasticMicroVelocityGradientdChi( dPlasticMicroLdChi );

                set_previousdPlasticMicroVelocityGradientdChin( dPlasticMicroLdChin );

                set_previousdPlasticMacroVelocityGradientdStateVariables( dPlasticMacroLdISVs );

                set_previousdPlasticMicroVelocityGradientdStateVariables( dPlasticMicroLdISVs );

                set_previousdPlasticGradientMicroVelocityGradientdMicroStress( dPlasticGradientMicroLdMicroStress );

                set_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( dPlasticGradientMicroLdHigherOrderStress );

                set_previousdPlasticGradientMicroVelocityGradientdF( dPlasticGradientMicroLdF );

                set_previousdPlasticGradientMicroVelocityGradientdFn( dPlasticGradientMicroLdFn );

                set_previousdPlasticGradientMicroVelocityGradientdChi( dPlasticGradientMicroLdChi );

                set_previousdPlasticGradientMicroVelocityGradientdChin( dPlasticGradientMicroLdChin );

                set_previousdPlasticGradientMicroVelocityGradientdGradChi( dPlasticGradientMicroLdGradChi );

                set_previousdPlasticGradientMicroVelocityGradientdGradChin( dPlasticGradientMicroLdGradChin );

                set_previousdPlasticGradientMicroVelocityGradientdStateVariables( dPlasticGradientMicroLdISVs );

            }
            else{

                set_plasticMacroVelocityGradient( macroVelocityGradient );

                set_plasticMicroVelocityGradient( microVelocityGradient );

                set_plasticGradientMicroVelocityGradient( gradientMicroVelocityGradient );

                set_dPlasticMacroVelocityGradientdMacroStress( dPlasticMacroLdMacroStress );

                set_dPlasticMacroVelocityGradientdMicroStress( dPlasticMacroLdMicroStress );

                set_dPlasticMicroVelocityGradientdMicroStress( dPlasticMicroLdMicroStress );

                set_dPlasticMacroVelocityGradientdF( dPlasticMacroLdF );

                set_dPlasticMacroVelocityGradientdFn( dPlasticMacroLdFn );

                set_dPlasticMicroVelocityGradientdF( dPlasticMicroLdF );

                set_dPlasticMicroVelocityGradientdFn( dPlasticMicroLdFn );

                set_dPlasticMicroVelocityGradientdChi( dPlasticMicroLdChi );

                set_dPlasticMicroVelocityGradientdChin( dPlasticMicroLdChin );

                set_dPlasticMacroVelocityGradientdStateVariables( dPlasticMacroLdISVs );

                set_dPlasticMicroVelocityGradientdStateVariables( dPlasticMicroLdISVs );

                set_dPlasticGradientMicroVelocityGradientdMicroStress( dPlasticGradientMicroLdMicroStress );

                set_dPlasticGradientMicroVelocityGradientdHigherOrderStress( dPlasticGradientMicroLdHigherOrderStress );

                set_dPlasticGradientMicroVelocityGradientdF( dPlasticGradientMicroLdF );

                set_dPlasticGradientMicroVelocityGradientdFn( dPlasticGradientMicroLdFn );

                set_dPlasticGradientMicroVelocityGradientdChi( dPlasticGradientMicroLdChi );

                set_dPlasticGradientMicroVelocityGradientdChin( dPlasticGradientMicroLdChin );

                set_dPlasticGradientMicroVelocityGradientdGradChi( dPlasticGradientMicroLdGradChi );

                set_dPlasticGradientMicroVelocityGradientdGradChin( dPlasticGradientMicroLdGradChin );

                set_dPlasticGradientMicroVelocityGradientdStateVariables( dPlasticGradientMicroLdISVs );

            }

        }

        void residual::setUpdatedPlasticDeformationGradient( ){
            /*!
             * Set the updated plastic deformation gradient
             */

            setPlasticDeformation( );

        }

        void residual::setUpdatedPlasticMicroDeformation( ){
            /*!
             * Set the updated plastic micro deformation
             */

            setPlasticDeformation( );

        }

        void residual::setUpdatedPlasticGradientMicroDeformation( ){
            /*!
             * Set the updated plastic gradient of the micro deformation
             */

            setPlasticDeformation( );

        }

        void residual::setPlasticDeformation( ){
            /*!
             * Set all of the plastic deformations
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            floatVector updatedPlasticDeformationGradient;

            floatVector updatedPlasticMicroDeformation;

            floatVector updatedPlasticGradientMicroDeformation;

            const floatVector previousPlasticDeformationGradient      = floatVector( hydra->get_previousConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector previousPlasticMicroDeformation         = floatVector( hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector previousPlasticGradientMicroDeformation = floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * ( plasticConfigurationIndex + 1 ) );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                evolvePlasticDeformation( *hydra->getDeltaTime( ),
                                          *get_plasticMacroVelocityGradient( ),
                                          *get_plasticMicroVelocityGradient( ),
                                          *get_plasticGradientMicroVelocityGradient( ),
                                          previousPlasticDeformationGradient,
                                          previousPlasticMicroDeformation,
                                          previousPlasticGradientMicroDeformation,
                                          *get_previousPlasticMacroVelocityGradient( ),
                                          *get_previousPlasticMicroVelocityGradient( ),
                                          *get_previousPlasticGradientMicroVelocityGradient( ),
                                          updatedPlasticDeformationGradient,
                                          updatedPlasticMicroDeformation,
                                          updatedPlasticGradientMicroDeformation,
                                          *getIntegrationParameter( ),
                                          *getIntegrationParameter( ),
                                          *getIntegrationParameter( ) );
            )

            set_updatedPlasticDeformationGradient( updatedPlasticDeformationGradient );

            set_updatedPlasticMicroDeformation( updatedPlasticMicroDeformation );

            set_updatedPlasticGradientMicroDeformation( updatedPlasticGradientMicroDeformation );

        }

        void residual::setdUpdatedPlasticDeformationGradientdMacroStress( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the macro stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticDeformationGradientdPreviousMacroStress( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the previous macro stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticDeformationGradientdMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the micro stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticDeformationGradientdPreviousMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the previous micro stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticDeformationGradientdF( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the deformation gradient
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticDeformationGradientdPreviousF( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the previous deformation gradient
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticDeformationGradientdFn( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the sub deformation gradients
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticDeformationGradientdPreviousFn( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the previous sub deformation gradients
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticDeformationGradientdStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the state variables
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticDeformationGradientdPreviousStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic deformation gradient w.r.t. the previous state variables
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the micro stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous micro stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdF( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the deformation gradient
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousF( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous deformation gradient
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdFn( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the sub deformation gradients
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousFn( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous sub deformation gradients
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdChi( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the micro deformation
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousChi( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous micro deformation
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdChin( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the sub micro deformations
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousChin( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous sub micro deformations
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticMicroDeformationdStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the state variables
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticMicroDeformationdPreviousStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic micro deformation w.r.t. the previous state variables
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdMacroStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the macro stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousMacroStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous macro stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the micro stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousMicroStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous micro stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdHigherOrderStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the higher order stress
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous higher order stress
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdF( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the deformation gradients
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousF( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous deformation gradients
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdFn( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the sub deformation gradients
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousFn( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous sub deformation gradients
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdChi( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the micro deformation
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousChi( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous micro deformation
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdChin( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the sub micro deformations
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousChin( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous sub micro deformations
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdGradChi( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the spatial gradient of the micro deformation
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousGradChi( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous spatial gradient of the micro deformation
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdGradChin( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the local spatial gradient of the sub micro deformations
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousGradChin( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous local spatial gradient of the sub micro deformations
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the state variables
             */

            setPlasticDeformationJacobians( false );

        }

        void residual::setdUpdatedPlasticGradientMicroDeformationdPreviousStateVariables( ){
            /*!
             * Set the jacobian of the updated plastic gradient micro deformation w.r.t. the previous state variables
             */

            setPlasticDeformationJacobians( true );

        }

        void residual::setPlasticDeformationJacobians( const bool addPreviousGradients  ){
            /*!
             * Set all of the plastic deformations
             *
              * \param addPreviousGradients: Flag for whether to compute the previous gradients
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            floatVector updatedPlasticDeformationGradient;

            floatVector updatedPlasticMicroDeformation;

            floatVector updatedPlasticGradientMicroDeformation;

            const floatVector previousPlasticDeformationGradient      = floatVector( hydra->get_previousConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector previousPlasticMicroDeformation         = floatVector( hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector previousPlasticGradientMicroDeformation = floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
                                                                                     hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * ( plasticConfigurationIndex + 1 ) );

            floatVector dPlasticFdPlasticMacroL;

            floatVector dPlasticMicroDeformationdPlasticMicroL;

            floatVector dPlasticGradientMicroDeformationdPlasticMacroL;

            floatVector dPlasticGradientMicroDeformationdPlasticMicroL;

            floatVector dPlasticGradientMicroDeformationdPlasticGradientMicroL;

            floatVector dPlasticFdPreviousPlasticF;

            floatVector dPlasticFdPreviousPlasticMacroL;

            if ( addPreviousGradients ){

                floatVector dPlasticFdPreviousPlasticF;
                floatVector dPlasticFdPreviousPlasticMacroL;
                floatVector dPlasticMicroDeformationdPreviousPlasticMicroDeformation;
                floatVector dPlasticMicroDeformationdPreviousPlasticMicroL;
                floatVector dPlasticGradientMicroDeformationdPreviousPlasticMicroDeformation;
                floatVector dPlasticGradientMicroDeformationdPreviousPlasticMicroGradient;
                floatVector dPlasticGradientMicroDeformationdPreviousPlasticMacroL;
                floatVector dPlasticGradientMicroDeformationdPreviousPlasticMicroL;
                floatVector dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL;

                TARDIGRADE_ERROR_TOOLS_CATCH(
                    evolvePlasticDeformation( *hydra->getDeltaTime( ),
                                              *get_plasticMacroVelocityGradient( ),
                                              *get_plasticMicroVelocityGradient( ),
                                              *get_plasticGradientMicroVelocityGradient( ),
                                              previousPlasticDeformationGradient,
                                              previousPlasticMicroDeformation,
                                              previousPlasticGradientMicroDeformation,
                                              *get_previousPlasticMacroVelocityGradient( ),
                                              *get_previousPlasticMicroVelocityGradient( ),
                                              *get_previousPlasticGradientMicroVelocityGradient( ),
                                              updatedPlasticDeformationGradient,
                                              updatedPlasticMicroDeformation,
                                              updatedPlasticGradientMicroDeformation,
                                              dPlasticFdPlasticMacroL,
                                              dPlasticMicroDeformationdPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPlasticMacroL,
                                              dPlasticGradientMicroDeformationdPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                              dPlasticFdPreviousPlasticF,
                                              dPlasticFdPreviousPlasticMacroL,
                                              dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                              dPlasticMicroDeformationdPreviousPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPreviousPlasticMicroDeformation,
                                              dPlasticGradientMicroDeformationdPreviousPlasticMicroGradient,
                                              dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                              dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                              *getIntegrationParameter( ),
                                              *getIntegrationParameter( ),
                                              *getIntegrationParameter( ) );
                )

                set_dUpdatedPlasticDeformationGradientdPreviousMacroStress( tardigradeVectorTools::matrixMultiply( dPlasticFdPreviousPlasticMacroL,
                                                                                                                   *get_previousdPlasticMacroVelocityGradientdMacroStress( ),
                                                                                                                   sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dUpdatedPlasticDeformationGradientdPreviousMicroStress( tardigradeVectorTools::matrixMultiply( dPlasticFdPreviousPlasticMacroL,
                                                                                                                   *get_previousdPlasticMacroVelocityGradientdMicroStress( ),
                                                                                                                   sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dUpdatedPlasticDeformationGradientdPreviousF( tardigradeVectorTools::matrixMultiply( dPlasticFdPreviousPlasticMacroL,
                                                                                                         *get_previousdPlasticMacroVelocityGradientdF( ),
                                                                                                         sot_dim, sot_dim, sot_dim, sot_dim ) );

                floatVector dUpdatedPlasticFdPreviousFn = tardigradeVectorTools::matrixMultiply( dPlasticFdPreviousPlasticMacroL,
                                                                                                 *get_previousdPlasticMacroVelocityGradientdFn( ),
                                                                                                 sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

                unsigned int offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dUpdatedPlasticFdPreviousFn[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticFdPreviousPlasticF[ sot_dim * i + j ];

                    }

                }

                set_dUpdatedPlasticDeformationGradientdPreviousFn( dUpdatedPlasticFdPreviousFn );

                set_dUpdatedPlasticDeformationGradientdPreviousStateVariables( tardigradeVectorTools::matrixMultiply( dPlasticFdPreviousPlasticMacroL,
                                                                                                                      *get_previousdPlasticMacroVelocityGradientdStateVariables( ),
                                                                                                                      sot_dim, sot_dim, sot_dim, num_isvs ) );

                set_dUpdatedPlasticMicroDeformationdPreviousMicroStress( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                                *get_previousdPlasticMicroVelocityGradientdMicroStress( ),
                                                                                                                sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dUpdatedPlasticMicroDeformationdPreviousF( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                      *get_previousdPlasticMicroVelocityGradientdF( ),
                                                                                                      sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dUpdatedPlasticMicroDeformationdPreviousFn( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                       *get_previousdPlasticMicroVelocityGradientdFn( ),
                                                                                                       sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dUpdatedPlasticMicroDeformationdPreviousChi( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                        *get_previousdPlasticMicroVelocityGradientdChi( ),
                                                                                                        sot_dim, sot_dim, sot_dim, sot_dim ) );

                floatVector dUpdatedPlasticMicroDeformationdPreviousChin = tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                                  *get_previousdPlasticMicroVelocityGradientdChin( ),
                                                                                                                  sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dUpdatedPlasticMicroDeformationdPreviousChin[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticMicroDeformationdPreviousPlasticMicroDeformation[ sot_dim * i + j ];

                    }

                }

                set_dUpdatedPlasticMicroDeformationdPreviousChin( dUpdatedPlasticMicroDeformationdPreviousChin );

                set_dUpdatedPlasticMicroDeformationdPreviousStateVariables( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPreviousPlasticMicroL,
                                                                                                                   *get_previousdPlasticMicroVelocityGradientdStateVariables( ),
                                                                                                                   sot_dim, sot_dim, sot_dim, num_isvs ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                                                                                                          *get_previousdPlasticMacroVelocityGradientdMacroStress( ),
                                                                                                                          tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                                                                                                          *get_previousdPlasticMacroVelocityGradientdMicroStress( ),
                                                                                                                          tot_dim, sot_dim, sot_dim, sot_dim )
                                                                                 + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                          *get_previousdPlasticMicroVelocityGradientdMicroStress( ),
                                                                                                                          tot_dim, sot_dim, sot_dim, sot_dim )
                                                                                 + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                          *get_previousdPlasticGradientMicroVelocityGradientdMicroStress( ),
                                                                                                                          tot_dim, tot_dim, tot_dim, sot_dim ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                              *get_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( ),
                                                                                                                              tot_dim, tot_dim, tot_dim, tot_dim ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousF(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                                                                                                *get_previousdPlasticMacroVelocityGradientdF( ),
                                                                                                                 tot_dim, sot_dim, sot_dim, sot_dim )
                                                                       + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                *get_previousdPlasticMicroVelocityGradientdF( ),
                                                                                                                tot_dim, sot_dim, sot_dim, sot_dim )
                                                                       + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                *get_previousdPlasticGradientMicroVelocityGradientdF( ),
                                                                                                                tot_dim, tot_dim, tot_dim, sot_dim ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousFn(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                                                                                                 *get_previousdPlasticMacroVelocityGradientdFn( ),
                                                                                                                 tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                        + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                 *get_previousdPlasticMicroVelocityGradientdFn( ),
                                                                                                                 tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                        + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                 *get_previousdPlasticGradientMicroVelocityGradientdFn( ),
                                                                                                                 tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousChi(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                  *get_previousdPlasticMicroVelocityGradientdChi( ),
                                                                                                                  tot_dim, sot_dim, sot_dim, sot_dim )
                                                                         + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                  *get_previousdPlasticGradientMicroVelocityGradientdChi( ),
                                                                                                                  tot_dim, tot_dim, tot_dim, sot_dim ) );

                floatVector dUpdatedPlasticGradientMicroDeformationdPreviousChin =  tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                           *get_previousdPlasticMicroVelocityGradientdChin( ),
                                                                                                                           tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                  + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                           *get_previousdPlasticGradientMicroVelocityGradientdChin( ),
                                                                                                                           tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < tot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dUpdatedPlasticGradientMicroDeformationdPreviousChin[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticGradientMicroDeformationdPreviousPlasticMicroDeformation[ sot_dim * i + j ];

                    }

                }

                set_dUpdatedPlasticGradientMicroDeformationdPreviousChin( dUpdatedPlasticGradientMicroDeformationdPreviousChin );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi( tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                    *get_previousdPlasticGradientMicroVelocityGradientdGradChi( ),
                                                                                                                    tot_dim, tot_dim, tot_dim, tot_dim ) );

                floatVector dUpdatedPlasticGradientMicroDeformationdPreviousGradChin = tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                              *get_previousdPlasticGradientMicroVelocityGradientdGradChin( ),
                                                                                                                              tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * tot_dim;

                for ( unsigned int i = 0; i < tot_dim; i++ ){

                    for ( unsigned int j = 0; j < tot_dim; j++ ){

                        dUpdatedPlasticGradientMicroDeformationdPreviousGradChin[ ( num_configs - 1 ) * tot_dim * i + j + offset ] += dPlasticGradientMicroDeformationdPreviousPlasticMicroGradient[ tot_dim * i  + j ];

                    }

                }

                set_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin( dUpdatedPlasticGradientMicroDeformationdPreviousGradChin );

                set_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMacroL,
                                                                                                                             *get_previousdPlasticMacroVelocityGradientdStateVariables( ),
                                                                                                                             tot_dim, sot_dim, sot_dim, num_isvs )
                                                                                    + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticMicroL,
                                                                                                                             *get_previousdPlasticMicroVelocityGradientdStateVariables( ),
                                                                                                                             tot_dim, sot_dim, sot_dim, num_isvs )
                                                                                    + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL,
                                                                                                                             *get_previousdPlasticGradientMicroVelocityGradientdStateVariables( ),
                                                                                                                             tot_dim, tot_dim, tot_dim, num_isvs ) );

            }
            else{
                TARDIGRADE_ERROR_TOOLS_CATCH(
                    evolvePlasticDeformation( *hydra->getDeltaTime( ),
                                              *get_plasticMacroVelocityGradient( ),
                                              *get_plasticMicroVelocityGradient( ),
                                              *get_plasticGradientMicroVelocityGradient( ),
                                              previousPlasticDeformationGradient,
                                              previousPlasticMicroDeformation,
                                              previousPlasticGradientMicroDeformation,
                                              *get_previousPlasticMacroVelocityGradient( ),
                                              *get_previousPlasticMicroVelocityGradient( ),
                                              *get_previousPlasticGradientMicroVelocityGradient( ),
                                              updatedPlasticDeformationGradient,
                                              updatedPlasticMicroDeformation,
                                              updatedPlasticGradientMicroDeformation,
                                              dPlasticFdPlasticMacroL,
                                              dPlasticMicroDeformationdPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPlasticMacroL,
                                              dPlasticGradientMicroDeformationdPlasticMicroL,
                                              dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                              *getIntegrationParameter( ),
                                              *getIntegrationParameter( ),
                                              *getIntegrationParameter( ) );
                )
            }

            set_updatedPlasticDeformationGradient( updatedPlasticDeformationGradient );

            set_updatedPlasticMicroDeformation( updatedPlasticMicroDeformation );

            set_updatedPlasticGradientMicroDeformation( updatedPlasticGradientMicroDeformation );

            set_dUpdatedPlasticDeformationGradientdMacroStress( tardigradeVectorTools::matrixMultiply( dPlasticFdPlasticMacroL,
                                                                                                       *get_dPlasticMacroVelocityGradientdMacroStress( ),
                                                                                                       sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticDeformationGradientdMicroStress( tardigradeVectorTools::matrixMultiply( dPlasticFdPlasticMacroL,
                                                                                                       *get_dPlasticMacroVelocityGradientdMicroStress( ),
                                                                                                       sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticDeformationGradientdF( tardigradeVectorTools::matrixMultiply( dPlasticFdPlasticMacroL,
                                                                                             *get_dPlasticMacroVelocityGradientdF( ),
                                                                                             sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticDeformationGradientdFn( tardigradeVectorTools::matrixMultiply( dPlasticFdPlasticMacroL,
                                                                                              *get_dPlasticMacroVelocityGradientdFn( ),
                                                                                              sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dUpdatedPlasticDeformationGradientdStateVariables( tardigradeVectorTools::matrixMultiply( dPlasticFdPlasticMacroL,
                                                                                                          *get_dPlasticMacroVelocityGradientdStateVariables( ),
                                                                                                          sot_dim, sot_dim, sot_dim, num_isvs ) );

            set_dUpdatedPlasticMicroDeformationdMicroStress( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                                    *get_dPlasticMicroVelocityGradientdMicroStress( ),
                                                                                                    sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticMicroDeformationdF( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                          *get_dPlasticMicroVelocityGradientdF( ),
                                                                                          sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticMicroDeformationdFn( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                           *get_dPlasticMicroVelocityGradientdFn( ),
                                                                                           sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dUpdatedPlasticMicroDeformationdChi( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                            *get_dPlasticMicroVelocityGradientdChi( ), sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticMicroDeformationdChin( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                             *get_dPlasticMicroVelocityGradientdChin( ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dUpdatedPlasticMicroDeformationdStateVariables( tardigradeVectorTools::matrixMultiply( dPlasticMicroDeformationdPlasticMicroL,
                                                                                                       *get_dPlasticMicroVelocityGradientdStateVariables( ),
                                                                                                       sot_dim, sot_dim, sot_dim, num_isvs ) );

            set_dUpdatedPlasticGradientMicroDeformationdMacroStress(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMacroL,
                                                                                                              *get_dPlasticMacroVelocityGradientdMacroStress( ),
                                                                                                              tot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdMicroStress(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMacroL,
                                                                                                              *get_dPlasticMacroVelocityGradientdMicroStress( ),
                                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                              *get_dPlasticMicroVelocityGradientdMicroStress( ),
                                                                                                              tot_dim, sot_dim, sot_dim, sot_dim )
                                                                     + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                              *get_dPlasticGradientMicroVelocityGradientdMicroStress( ),
                                                                                                              tot_dim, tot_dim, tot_dim, sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress( tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                                  *get_dPlasticGradientMicroVelocityGradientdHigherOrderStress( ),
                                                                                                                  tot_dim, tot_dim, tot_dim, tot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdF(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMacroL,
                                                                                                    *get_dPlasticMacroVelocityGradientdF( ),
                                                                                                    tot_dim, sot_dim, sot_dim, sot_dim )
                                                           + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                    *get_dPlasticMicroVelocityGradientdF( ),
                                                                                                    tot_dim, sot_dim, sot_dim, sot_dim )
                                                           + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                    *get_dPlasticGradientMicroVelocityGradientdF( ),
                                                                                                    tot_dim, tot_dim, tot_dim, sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdFn(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMacroL,
                                                                                                     *get_dPlasticMacroVelocityGradientdFn( ),
                                                                                                     tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                            + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                     *get_dPlasticMicroVelocityGradientdFn( ),
                                                                                                     tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                            + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                     *get_dPlasticGradientMicroVelocityGradientdFn( ),
                                                                                                     tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdChi(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                      *get_dPlasticMicroVelocityGradientdChi( ),
                                                                                                      tot_dim, sot_dim, sot_dim, sot_dim )
                                                             + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                      *get_dPlasticGradientMicroVelocityGradientdChi( ),
                                                                                                      tot_dim, tot_dim, tot_dim, sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdChin(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                       *get_dPlasticMicroVelocityGradientdChin( ),
                                                                                                       tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                              + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                       *get_dPlasticGradientMicroVelocityGradientdChin( ),
                                                                                                       tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdGradChi( tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                        *get_dPlasticGradientMicroVelocityGradientdGradChi( ),
                                                                                                        tot_dim, tot_dim, tot_dim, tot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdGradChin( tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                         *get_dPlasticGradientMicroVelocityGradientdGradChin( ),
                                                                                                         tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ) );

            set_dUpdatedPlasticGradientMicroDeformationdStateVariables(   tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMacroL,
                                                                                                                 *get_dPlasticMacroVelocityGradientdStateVariables( ),
                                                                                                                 tot_dim, sot_dim, sot_dim, num_isvs )
                                                                        + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticMicroL,
                                                                                                                 *get_dPlasticMicroVelocityGradientdStateVariables( ),
                                                                                                                 tot_dim, sot_dim, sot_dim, num_isvs )
                                                                        + tardigradeVectorTools::matrixMultiply( dPlasticGradientMicroDeformationdPlasticGradientMicroL,
                                                                                                                 *get_dPlasticGradientMicroVelocityGradientdStateVariables( ),
                                                                                                                 tot_dim, tot_dim, tot_dim, num_isvs ) );

        }

        void residual::setStateVariableResiduals( ){
            /*!
             * Set the state variable residuals
             *
             * We define these residuals as
             *
             * \f$R = \left\langle f \right\rangle - \dot{\gamma} \left\langle -f \right\rangle\f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const floatVector *updatedPlasticStrainLikeISVs = get_updatedPlasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const floatVector *microGradientYield = get_microGradientYield( );

            floatVector residual( get_plasticStateVariables( )->size( ), 0 );

            floatType macroMac  = tardigradeConstitutiveTools::mac( *macroYield );
            tardigradeConstitutiveTools::mac( -( *macroYield ) );

            floatType microMac  = tardigradeConstitutiveTools::mac( *microYield );
            tardigradeConstitutiveTools::mac( -( *microYield ) );

            floatVector microGradientMac( microGradientYield->size( ), 0 );
            floatVector nMicroGradientMac( microGradientYield->size( ), 0 );

            floatType macNegMacroGamma;
            floatType macNegMicroGamma;

            floatVector macNegMicroGradientGamma( 3, 0 );

            if ( *useWeakenedMacaulay( ) ){

                macroMac  = weakMac( *macroYield, *getWeakenedMacaulayParameter( ) );
                weakMac( -( *macroYield ), *getWeakenedMacaulayParameter( ) );

                microMac  = weakMac( *microYield, *getWeakenedMacaulayParameter( ) );
                weakMac( -( *microYield ), *getWeakenedMacaulayParameter( ) );

                macNegMacroGamma = weakMac( -( *plasticMultipliers )[ 0 ] , *getWeakenedMacaulayParameter( ) );
                macNegMicroGamma = weakMac( -( *plasticMultipliers )[ 1 ] , *getWeakenedMacaulayParameter( ) );

                for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){
                    unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );
                    microGradientMac[ index ]  = weakMac( *y, *getWeakenedMacaulayParameter( ) );
                    nMicroGradientMac[ index ] = weakMac( -( *y ), *getWeakenedMacaulayParameter( ) );
                    macNegMicroGradientGamma[ index ] = weakMac( -( *plasticMultipliers )[ index + 2 ], *getWeakenedMacaulayParameter( ) );
                }
            }
            else{

                macroMac  = tardigradeConstitutiveTools::mac( *macroYield );
                tardigradeConstitutiveTools::mac( -( *macroYield ) );

                microMac  = tardigradeConstitutiveTools::mac( *microYield );
                tardigradeConstitutiveTools::mac( -( *microYield ) );

                macNegMacroGamma = tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ 0 ] );
                macNegMicroGamma = tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ 1 ] );

                for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){
                    unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );
                    microGradientMac[ index ] = tardigradeConstitutiveTools::mac( *y );
                    nMicroGradientMac[ index ] = tardigradeConstitutiveTools::mac( -( *y ) );
                    macNegMicroGradientGamma[ index ] = tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ index + 2 ] );
                }

            }

            // Set the terms associated with the yield surface
            residual[ 0 ] = macroMac + ( *plasticMultipliers )[ 0 ] * ( *macroYield ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMacroGamma;

            residual[ 1 ] = microMac + ( *plasticMultipliers )[ 1 ] * ( *microYield ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMicroGamma;

            for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){

                unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );

                residual[ index + 2 ]
                    = microGradientMac[ index ] + ( *plasticMultipliers )[ index + 2 ] * ( *y ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMicroGradientGamma[ index ];

            }

            // Set the terms associated with the strain-like ISV evolution
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                residual[ numPlasticMultipliers + i ] = ( *updatedPlasticStrainLikeISVs )[ i ] - ( *plasticStrainLikeISVs )[ i ];

            }

            set_stateVariableResiduals( residual );

        }

        void residual::setStateVariableJacobians( ){
            /*!
             * Set the state variable residual jacobians
             *
             * We define these residuals as
             *
             * \f$R = \left\langle f \right\rangle - \dot{\gamma} \left\langle -f \right\rangle\f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const unsigned int numUnknowns = hydra->getUnknownVector( )->size( );

            const unsigned int numISVs = get_plasticStateVariables( )->size( );

            const floatVector *dMacroYielddStress                 = get_dMacroYielddStress( );

            const floatVector *dMacroYielddFn                     = get_dMacroYielddFn( );

            const floatVector *dMacroYielddStateVariables         = get_dMacroYielddStateVariables( );

            const floatVector *dMicroYielddStress                 = get_dMicroYielddStress( );

            const floatVector *dMicroYielddFn                     = get_dMicroYielddFn( );

            const floatVector *dMicroYielddStateVariables         = get_dMicroYielddStateVariables( );

            const floatVector *dMicroGradientYielddStress         = get_dMicroGradientYielddStress( );

            const floatVector *dMicroGradientYielddFn             = get_dMicroGradientYielddFn( );

            const floatVector *dMicroGradientYielddChin           = get_dMicroGradientYielddChin( );

            const floatVector *dMicroGradientYielddStateVariables = get_dMicroGradientYielddStateVariables( );

            const floatVector *dUpdatedPlasticStrainLikeISVsdStateVariables = get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const floatVector *microGradientYield = get_microGradientYield( );

            floatVector jacobian( numISVs * numUnknowns, 0 );

            // Stress Jacobians
            floatType dMacroMacdx, dMicroMacdx;
            floatType ndMacroMacdx, ndMicroMacdx;

            floatVector microGradientMac( numPlasticMultipliers - 2 );
            floatVector dMicroGradientMacdx( numPlasticMultipliers - 2 );

            floatVector nMicroGradientMac( numPlasticMultipliers - 2 );
            floatVector ndMicroGradientMacdx( numPlasticMultipliers - 2 );

            floatType dMacNegMacroGammadGamma;
            floatType dMacNegMicroGammadGamma;

            floatVector dMacNegMicroGradientGammadGamma( 3, 0 );

            if ( *useWeakenedMacaulay( ) ){

                weakMac( *macroYield, *getWeakenedMacaulayParameter( ), dMacroMacdx );

                weakMac( *microYield, *getWeakenedMacaulayParameter( ), dMicroMacdx );

                weakMac( -( *macroYield ), *getWeakenedMacaulayParameter( ), ndMacroMacdx );

                weakMac( -( *microYield ), *getWeakenedMacaulayParameter( ), ndMicroMacdx );

                weakMac( -( *plasticMultipliers )[ 0 ], *getWeakenedMacaulayParameter( ), dMacNegMacroGammadGamma );

                weakMac( -( *plasticMultipliers )[ 1 ], *getWeakenedMacaulayParameter( ), dMacNegMicroGammadGamma );

                for ( unsigned int i = 0; i < ( numPlasticMultipliers - 2 ); i++ ){

                    microGradientMac[ i ]  = weakMac(  ( *microGradientYield )[ i ], *getWeakenedMacaulayParameter( ), dMicroGradientMacdx[ i ] );

                    nMicroGradientMac[ i ] = weakMac( -( *microGradientYield )[ i ], *getWeakenedMacaulayParameter( ), ndMicroGradientMacdx[ i ] );

                    weakMac( -( *plasticMultipliers )[ i + 2 ], *getWeakenedMacaulayParameter( ), dMacNegMicroGradientGammadGamma[ i ] );

                }

            }
            else{

                tardigradeConstitutiveTools::mac( *macroYield, dMacroMacdx );

                tardigradeConstitutiveTools::mac( *microYield, dMicroMacdx );

                tardigradeConstitutiveTools::mac( -( *macroYield ), ndMacroMacdx );

                tardigradeConstitutiveTools::mac( -( *microYield ), ndMicroMacdx );

                tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ 0 ], dMacNegMacroGammadGamma );

                tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ 1 ], dMacNegMicroGammadGamma );

                for ( unsigned int i = 0; i < ( numPlasticMultipliers - 2 ); i++ ){

                    microGradientMac[ i ]  = tardigradeConstitutiveTools::mac(  ( *microGradientYield )[ i ],  dMicroGradientMacdx[ i ] );

                    nMicroGradientMac[ i ] = tardigradeConstitutiveTools::mac( -( *microGradientYield )[ i ], ndMicroGradientMacdx[ i ] );

                    tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ i + 2 ], dMacNegMicroGradientGammadGamma[ i ] );

                }

            }

            unsigned int offset = numSecondOrderTensor;
            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                jacobian[ numUnknowns * 0 + j ] = ( dMacroMacdx  + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddStress )[ j ];

                jacobian[ numUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddStress )[ j ];

            }

            offset = 2 * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    jacobian[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddStress )[ numThirdOrderTensor * i + j ];

                }

            }

            // Sub-Deformation gradient jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor;
            for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                jacobian[ numUnknowns * 0 + j + offset ] = ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddFn )[ j ];

                jacobian[ numUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddFn )[ j ];

            } 

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }

            // Sub-Micro deformation jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor + ( numConfigurations - 1 ) * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }


            // State Variable Jacobians
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );

            jacobian[ numUnknowns * 0 + offset + 0 ] += *macroYield - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMacroGammadGamma;

            jacobian[ numUnknowns * 1 + offset + 1 ] += *microYield - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMicroGammadGamma;

            for ( unsigned int i = 0; i < dim; i++ ){

                jacobian[ numUnknowns * ( i + 2 ) + offset + i + 2 ] += ( *microGradientYield )[ i ] - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMicroGradientGammadGamma[ i ];

            }

            for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                jacobian[ numUnknowns * 0 + j + offset ] += ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddStateVariables )[ j ];

                jacobian[ numUnknowns * 1 + j + offset ] += ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddStateVariables )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                    jacobian[ numUnknowns * ( i + 2 ) + j + offset ] += ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddStateVariables )[ numISVs * i + j ];

                }

            }
            
            unsigned int row0 = numPlasticMultipliers;
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                jacobian[ numUnknowns * ( i + row0 ) + i + offset + numPlasticMultipliers ] -= 1;

                for ( auto j = getStateVariableIndices( )->begin( ); j != getStateVariableIndices( )->end( ); j++ ){

                    jacobian[ numUnknowns * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numISVs * i + ( unsigned int )( j - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            set_stateVariableJacobians( jacobian );

        }

        void residual::setdStateVariableResidualsdD( ){
            /*!
             * Set the state variable residuals derivatives w.r.t. the deformation measures
             *
             * We define these residuals as
             *
             * \f$R = \left\langle f \right\rangle - \dot{\gamma} \left\langle -f \right\rangle\f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int dim = hydra->getDimension( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            const floatVector *dMacroYielddF                      = get_dMacroYielddF( );

            const floatVector *dMicroYielddF                      = get_dMicroYielddF( );

            const floatVector *dMicroGradientYielddF              = get_dMicroGradientYielddF( );

            const floatVector *dMicroGradientYielddChi            = get_dMicroGradientYielddChi( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const floatVector *microGradientYield = get_microGradientYield( );

            floatVector dRdD( get_plasticStateVariables( )->size( ) * numConfigurationUnknowns, 0 );

            // Deformation gradient jacobians
            floatType dMacroMacdx, dMicroMacdx;
            floatType ndMacroMacdx, ndMicroMacdx;

            floatVector dMicroGradientMacdx( numPlasticMultipliers - 2 );

            floatVector ndMicroGradientMacdx( numPlasticMultipliers - 2 );

            if ( *useWeakenedMacaulay( ) ){

                weakMac( *macroYield, *getWeakenedMacaulayParameter( ), dMacroMacdx );

                weakMac( *microYield, *getWeakenedMacaulayParameter( ), dMicroMacdx );

                weakMac( -( *macroYield ), *getWeakenedMacaulayParameter( ), ndMacroMacdx );

                weakMac( -( *microYield ), *getWeakenedMacaulayParameter( ), ndMicroMacdx );

                for ( unsigned int i = 0; i < ( numPlasticMultipliers - 2 ); i++ ){

                    weakMac(  ( *microGradientYield )[ i ], *getWeakenedMacaulayParameter( ), dMicroGradientMacdx[ i ] );

                    weakMac( -( *microGradientYield )[ i ], *getWeakenedMacaulayParameter( ),ndMicroGradientMacdx[ i ] );

                }

            }
            else{

                tardigradeConstitutiveTools::mac( *macroYield, dMacroMacdx );

                tardigradeConstitutiveTools::mac( *microYield, dMicroMacdx );

                tardigradeConstitutiveTools::mac( -( *macroYield ), ndMacroMacdx );

                tardigradeConstitutiveTools::mac( -( *microYield ), ndMicroMacdx );

                for ( unsigned int i = 0; i < ( numPlasticMultipliers - 2 ); i++ ){

                    tardigradeConstitutiveTools::mac(  ( *microGradientYield )[ i ],  dMicroGradientMacdx[ i ] );

                    tardigradeConstitutiveTools::mac( -( *microGradientYield )[ i ], ndMicroGradientMacdx[ i ] );

                }

            }

            unsigned int offset = 0;
            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                dRdD[ numConfigurationUnknowns * 0 + j + offset ] = ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddF )[ j ];

                dRdD[ numConfigurationUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddF )[ j ];

            } 

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    dRdD[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddF )[ numSecondOrderTensor * i + j ];

                }

            }

            // Micro deformation jacobians
            offset = numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    dRdD[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddChi )[ numSecondOrderTensor * i + j ];

                }

            }

            set_dStateVariableResidualsdD( dRdD );

        }

        void residual::setdStateVariableResidualsdPreviousISVs( ){
            /*!
             * Set the derivatives of the state variable residuals with respect to the previous ISV vector
             *
             * We define these residuals as
             *
             * \f$R = \left\langle f \right\rangle - \dot{\gamma} \left\langle -f \right\rangle\f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const floatVector *dUpdatedPlasticStrainLikeISVsdStateVariables = get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const unsigned int numPlasticISVs = get_plasticStateVariables( )->size( );

            const unsigned int numISVs = hydra->getPreviousStateVariables( )->size( );

            floatVector dRdPreviousISVs( numPlasticISVs * numISVs, 0 );

            // Stress Jacobians
            unsigned int row0 = numPlasticMultipliers;
            unsigned int offset = ( numConfigurations - 1 ) * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            std::vector< unsigned int > stateVariableIndices = *getStateVariableIndices( );

            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                dRdPreviousISVs[ numISVs * ( i + row0 ) + stateVariableIndices[ i ] + offset + numPlasticMultipliers ] += 1;

                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    dRdPreviousISVs[ numISVs * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numPlasticISVs * i + ( unsigned int )( j - stateVariableIndices.begin( ) ) ];

                }

            }

            set_dStateVariableResidualsdPreviousISVs( dRdPreviousISVs );

        }

        void residual::setResidual( ){
            /*!
             * Set the residual equation
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );
 
            const unsigned int tot_dim = hydra->getTOTDimension( ); 

            const floatVector *updatedPlasticDeformationGradient;

            const floatVector *updatedPlasticMicroDeformation;

            const floatVector *updatedPlasticGradientMicroDeformation;

            const floatVector *stateVariableResiduals;

            // Get the trial plastic deformation measures
            unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            const floatVector plasticDeformationGradient      = floatVector( hydra->get_configurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                             hydra->get_configurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector plasticMicroDeformation         = floatVector( hydra->get_microConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                             hydra->get_microConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const floatVector plasticGradientMicroDeformation = floatVector( hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
                                                                             hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim * ( plasticConfigurationIndex + 1 ) );

            // Get the updated plastic deformation measures
            TARDIGRADE_ERROR_TOOLS_CATCH(
                updatedPlasticDeformationGradient = get_updatedPlasticDeformationGradient( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                updatedPlasticMicroDeformation = get_updatedPlasticMicroDeformation( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                updatedPlasticGradientMicroDeformation = get_updatedPlasticGradientMicroDeformation( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                stateVariableResiduals = get_stateVariableResiduals( );
            )

            floatVector residual = tardigradeVectorTools::appendVectors( { *updatedPlasticDeformationGradient      - plasticDeformationGradient,
                                                                           *updatedPlasticMicroDeformation         - plasticMicroDeformation,
                                                                           *updatedPlasticGradientMicroDeformation - plasticGradientMicroDeformation,
                                                                           *stateVariableResiduals } );

            setResidual( residual );

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian
             */

            const unsigned int numEquations = *getNumEquations( );

            const unsigned int numUnknowns  = hydra->getUnknownVector( )->size( );

            const unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            const std::vector< unsigned int > stateVariableIndices = *getStateVariableIndices( );

            const unsigned int numISVs = stateVariableIndices.size( );

            const floatVector *dUpdatedPlasticDeformationGradientdMacroStress;

            const floatVector *dUpdatedPlasticDeformationGradientdMicroStress;

            const floatVector *dUpdatedPlasticDeformationGradientdFn;

            const floatVector *dUpdatedPlasticDeformationGradientdStateVariables;

            const floatVector *dUpdatedPlasticMicroDeformationdMicroStress;

            const floatVector *dUpdatedPlasticMicroDeformationdFn;

            const floatVector *dUpdatedPlasticMicroDeformationdChin;

            const floatVector *dUpdatedPlasticMicroDeformationdStateVariables;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdMacroStress;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdMicroStress;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdHigherOrderStress;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdFn;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdChin;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdGradChin;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdStateVariables;

            const floatVector *stateVariableJacobians;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticDeformationGradientdMacroStress = get_dUpdatedPlasticDeformationGradientdMacroStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticDeformationGradientdMicroStress = get_dUpdatedPlasticDeformationGradientdMicroStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticDeformationGradientdFn = get_dUpdatedPlasticDeformationGradientdFn( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticDeformationGradientdStateVariables = get_dUpdatedPlasticDeformationGradientdStateVariables( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdMicroStress = get_dUpdatedPlasticMicroDeformationdMicroStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdFn = get_dUpdatedPlasticMicroDeformationdFn( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdChin = get_dUpdatedPlasticMicroDeformationdChin( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdStateVariables = get_dUpdatedPlasticMicroDeformationdStateVariables( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdMacroStress = get_dUpdatedPlasticGradientMicroDeformationdMacroStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdMicroStress = get_dUpdatedPlasticGradientMicroDeformationdMicroStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdHigherOrderStress = get_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdFn = get_dUpdatedPlasticGradientMicroDeformationdFn( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdChin = get_dUpdatedPlasticGradientMicroDeformationdChin( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdGradChin = get_dUpdatedPlasticGradientMicroDeformationdGradChin( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdStateVariables = get_dUpdatedPlasticGradientMicroDeformationdStateVariables( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                stateVariableJacobians = get_stateVariableJacobians( );
            )

            floatMatrix jacobian( numEquations, floatVector( numUnknowns, 0 ) );

            // Set the Jacobians of the second order plastic deformation measures
            for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                // Jacobians with respect to the trial stresses
                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    jacobian[ i                        ][ j                        ] += ( *dUpdatedPlasticDeformationGradientdMacroStress )[ numSecondOrderTensor * i + j ];

                    jacobian[ i                        ][ j + numSecondOrderTensor ] += ( *dUpdatedPlasticDeformationGradientdMicroStress )[ numSecondOrderTensor * i + j ];

                    jacobian[ i + numSecondOrderTensor ][ j + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdMicroStress )[ numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-deformation gradients
                jacobian[ i                        ][ i + 2 * numSecondOrderTensor + numThirdOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ i                        ][ j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticDeformationGradientdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                    jacobian[ i + numSecondOrderTensor ][ j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-micro deformation

                jacobian[ i + numSecondOrderTensor ][ i + numConfigurationUnknowns + numSecondOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ i + numSecondOrderTensor ][ j + numConfigurationUnknowns + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the state variables
                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    unsigned int col = ( unsigned int )( j - stateVariableIndices.begin( ) );

                    jacobian[ i                        ][ *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticDeformationGradientdStateVariables )[ numISVs * i + col ];

                    jacobian[ i + numSecondOrderTensor ][ *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticMicroDeformationdStateVariables )[ numISVs * i + col ];

                }

            }

            // Set the Jacobians of the third order plastic deformation measures
            for ( unsigned int i = 0; i < numThirdOrderTensor; i++ ){

                // Set the Jacobians with respect to the trial stresses
                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    jacobian[ i + 2 * numSecondOrderTensor ][ j                        ] = ( *dUpdatedPlasticGradientMicroDeformationdMacroStress )[ numSecondOrderTensor * i + j ];

                    jacobian[ i + 2 * numSecondOrderTensor ][ j + numSecondOrderTensor ] = ( *dUpdatedPlasticGradientMicroDeformationdMicroStress )[ numSecondOrderTensor * i + j ];

                }

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    jacobian[ i + 2 * numSecondOrderTensor ][ j + 2 * numSecondOrderTensor ] = ( *dUpdatedPlasticGradientMicroDeformationdHigherOrderStress )[ numThirdOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ i + 2 * numSecondOrderTensor ][ j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-micro deformations
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    jacobian[ i + 2 * numSecondOrderTensor ][ j + numConfigurationUnknowns + numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Set the jacobians with respect to the local spatial gradients of the sub-micro deformation
                jacobian[ i + 2 * numSecondOrderTensor ][ i + numConfigurationUnknowns + 2 * numSecondOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numThirdOrderTensor; j++ ){

                    jacobian[ i + 2 * numSecondOrderTensor ][ j + numConfigurationUnknowns + 2 * numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdGradChin )[ ( numConfigurations - 1 ) * numThirdOrderTensor * i + j ];

                }

                // Jacobians with respect to the state variables
                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    unsigned int col = ( unsigned int )( j - stateVariableIndices.begin( ) );

                    jacobian[ i + 2 * numSecondOrderTensor ][ *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticGradientMicroDeformationdStateVariables )[ numISVs * i + col ];

                }

            }

            // Set the Jacobians of the state variables
            for ( unsigned int row = 0; row < numISVs; row++ ){

                jacobian[ numConfigurationUnknowns + row ] = floatVector( stateVariableJacobians->begin( ) + numUnknowns * row,
                                                                          stateVariableJacobians->begin( ) + numUnknowns * ( row + 1 ) );

            }

            setJacobian( jacobian );

        }

        void residual::setdRdD( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation
             */

            const unsigned int numEquations = *getNumEquations( );

            const unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            const std::vector< unsigned int > stateVariableIndices = *getStateVariableIndices( );

            const unsigned int numISVs = stateVariableIndices.size( );

            const floatVector *dUpdatedPlasticDeformationGradientdF;

            const floatVector *dUpdatedPlasticMicroDeformationdF;

            const floatVector *dUpdatedPlasticMicroDeformationdChi;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdF;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdChi;

            const floatVector *dUpdatedPlasticGradientMicroDeformationdGradChi;

            const floatVector *dStateVariableResidualsdD;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticDeformationGradientdF = get_dUpdatedPlasticDeformationGradientdF( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdF = get_dUpdatedPlasticMicroDeformationdF( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticMicroDeformationdChi = get_dUpdatedPlasticMicroDeformationdChi( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdF = get_dUpdatedPlasticGradientMicroDeformationdF( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdChi = get_dUpdatedPlasticGradientMicroDeformationdChi( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dUpdatedPlasticGradientMicroDeformationdGradChi = get_dUpdatedPlasticGradientMicroDeformationdGradChi( );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                dStateVariableResidualsdD = get_dStateVariableResidualsdD( );
            )

            floatMatrix dRdD( numEquations, floatVector( *hydra->getConfigurationUnknownCount( ), 0 ) );

            // Set the Jacobians of the second order plastic deformation measures
            for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                // Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    dRdD[ i                        ][ j ] += ( *dUpdatedPlasticDeformationGradientdF )[ numSecondOrderTensor * i + j ];

                    dRdD[ i + numSecondOrderTensor ][ j ] += ( *dUpdatedPlasticMicroDeformationdF )[ numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-micro deformation
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    dRdD[ i + numSecondOrderTensor ][ j + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdChi )[ numSecondOrderTensor * i + j ];
 
                }

            }

            // Set the Jacobians of the third order plastic deformation measures
            for ( unsigned int i = 0; i < numThirdOrderTensor; i++ ){

                // Set the Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    dRdD[ i + 2 * numSecondOrderTensor ][ j ] += ( *dUpdatedPlasticGradientMicroDeformationdF )[ numSecondOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-micro deformations
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    dRdD[ i + 2 * numSecondOrderTensor ][ j + numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdChi )[ numSecondOrderTensor * i + j ];

                }

                // Set the jacobians with respect to the local spatial gradients of the sub-micro deformation
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numThirdOrderTensor; j++ ){

                    dRdD[ i + 2 * numSecondOrderTensor ][ j + 2 * numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdGradChi )[ numThirdOrderTensor * i + j ];

                }

            }

            // Set the Jacobians of the state variables
            for ( unsigned int row = 0; row < numISVs; row++ ){

                dRdD[ numConfigurationUnknowns + row ] = floatVector( dStateVariableResidualsdD->begin( ) + numConfigurationUnknowns * row,
                                                                      dStateVariableResidualsdD->begin( ) + numConfigurationUnknowns * ( row + 1 ) );

            }

            setdRdD( dRdD );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            setdRdT( floatVector( *getNumEquations( ), 0 ) );

        }

    }

}
