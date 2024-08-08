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
            TARDIGRADE_ERROR_TOOLS_CHECK( ( 0 <= frictionAngle ) && ( frictionAngle <= 1.570796 ), "The friction angle must be betwen 0 and pi / 2 not " + std::to_string( frictionAngle ) );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( abs( beta ) <= 1, "Beta must be between -1 and 1 not " + std::to_string( beta ) );
    
            //Compute the parameters
            parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );
    
            A = betaAngle * std::cos( frictionAngle );
    
            B = betaAngle * std::sin( frictionAngle );
    
        }

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &stressMeasure, const variableType &cohesion,
                                                               const secondOrderTensor &precedingDeformationGradient,
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
            secondOrderTensor rightCauchyGreen;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen ) );

            //Compute the decomposition of the stress
            variableType pressure;
            secondOrderTensor deviatoricReferenceStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( stressMeasure,
                                                                                                                                    rightCauchyGreen, deviatoricReferenceStress, pressure ) );

            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        }

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &referenceStressMeasure, const variableType &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, secondOrderTensor &dFdStress, variableType &dFdc,
                                                           secondOrderTensor &dFdPrecedingF, double tol ){
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
            secondOrderTensor rightCauchyGreen;
            fourthOrderTensor dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );
 
            //Compute the decomposition of the stress
            variableType pressure;
            secondOrderTensor deviatoricReferenceStress;
    
            fourthOrderTensor dDevStressdStress, dDevStressdRCG;
            secondOrderTensor dPressuredStress, dPressuredRCG;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );
   
            fourthOrderTensor dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim );
            secondOrderTensor dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG, dRCGdPrecedingF, 1, sot_dim, sot_dim, sot_dim );
 
            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );
    
            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Evaluate the jacobians
            secondOrderTensor devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );
    
            dFdStress  = BAngle * dPressuredStress;
            dFdStress += tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdStress, 1, sot_dim, sot_dim, sot_dim );
    
            dFdc = - AAngle;
    
            dFdPrecedingF  = BAngle * dPressuredPrecedingF;
            dFdPrecedingF += tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdPrecedingF, 1, sot_dim, sot_dim, sot_dim );
        }

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &stressMeasure, const variableType &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, secondOrderTensor &dFdStress, variableType &dFdc,
                                                           secondOrderTensor &dFdPrecedingF, fourthOrderTensor &d2FdStress2,
                                                           fourthOrderTensor &d2FdStressdPrecedingF, double tol ){
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
            constexpr unsigned int fot_dim = sot_dim * sot_dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH(  computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            secondOrderTensor rightCauchyGreen;
            fourthOrderTensor dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableType pressure;
            secondOrderTensor deviatoricReferenceStress;

            fourthOrderTensor dDevStressdStress, dDevStressdRCG;
            secondOrderTensor dPressuredStress, dPressuredRCG;

            sixthOrderTensor d2DevStressdStressdRCG;
            fourthOrderTensor d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            fourthOrderTensor dDevStressdPrecedingF = tardigradeVectorTools::matrixMultiply( dDevStressdRCG, dRCGdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim );
            secondOrderTensor dPressuredPrecedingF  = tardigradeVectorTools::matrixMultiply( dPressuredRCG, dRCGdPrecedingF, 1, sot_dim, sot_dim, sot_dim );

            sixthOrderTensor d2DevStressdStressdPrecedingF( sot_dim * sot_dim * sot_dim, 0 );

            fourthOrderTensor d2PressuredStressdPrecedingF( sot_dim * sot_dim, 0 );

            for ( unsigned int I = 0; I < deviatoricReferenceStress.size( ); I++ ){

                for ( unsigned int J = 0; J < deviatoricReferenceStress.size( ); J++ ){

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < rightCauchyGreen.size( ); L++ ){

                            d2DevStressdStressdPrecedingF[ sot_dim * sot_dim * I + sot_dim * J + L ]
                                += d2DevStressdStressdRCG[ sot_dim * sot_dim * I + sot_dim * J + K ] * dRCGdPrecedingF[ sot_dim * K + L ];

                        }

                        d2PressuredStressdPrecedingF[ sot_dim * I + K ]
                            += d2PressuredStressdRCG[ sot_dim * I + J ] * dRCGdPrecedingF[ sot_dim * J + K ];

                    }

                }

            }

            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Evaluate the jacobians
            secondOrderTensor devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

            dFdStress = BAngle * dPressuredStress;
            dFdStress += tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdStress, 1, sot_dim, sot_dim, sot_dim );

            dFdc = - AAngle;

            dFdPrecedingF = BAngle * dPressuredPrecedingF;
            dFdPrecedingF += tardigradeVectorTools::matrixMultiply( devStressDirection, dDevStressdPrecedingF, 1, sot_dim, sot_dim, sot_dim );

            //Evaluate the second-order jacobians
            fourthOrderTensor dDevStressDirectiondDevStress( sot_dim * sot_dim, 0 );
            for ( unsigned int i = 0; i < sot_dim; i++ ){ dDevStressDirectiondDevStress[ sot_dim * i + i ] = 1. / ( normDevStress + tol ); }
            dDevStressDirectiondDevStress -= tardigradeVectorTools::matrixMultiply( devStressDirection, devStressDirection, sot_dim, 1, 1, sot_dim ) / ( normDevStress + tol );

            d2FdStress2 = tardigradeVectorTools::matrixMultiply( dDevStressdStress, tardigradeVectorTools::matrixMultiply( dDevStressDirectiondDevStress, dDevStressdStress, sot_dim, sot_dim, sot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim, true, false );

            d2FdStressdPrecedingF  = BAngle * d2PressuredStressdPrecedingF;
            d2FdStressdPrecedingF += tardigradeVectorTools::matrixMultiply( tardigradeVectorTools::matrixMultiply( dDevStressDirectiondDevStress, dDevStressdStress, sot_dim, sot_dim, sot_dim, sot_dim ), dDevStressdPrecedingF, sot_dim, sot_dim, sot_dim, sot_dim, true, false );

            for ( unsigned int AB = 0; AB < sot_dim; AB++ ){
                for ( unsigned int IJKL = 0; IJKL < fot_dim; IJKL++ ){
                    d2FdStressdPrecedingF[ IJKL ] += devStressDirection[ AB ] * d2DevStressdStressdPrecedingF[ fot_dim * AB + IJKL ];
                }
            }

        }

        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor &stressMeasure,
                                                               const dimVector &cohesion,
                                                               const secondOrderTensor &precedingDeformationGradient,
                                                               const parameterType &frictionAngle, const parameterType &beta,
                                                               dimVector &yieldValue ){
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
            secondOrderTensor rightCauchyGreen;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen ) );

            //Compute the decomposition of the stress
            dimVector pressure;
            thirdOrderTensor deviatoricReferenceStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                 rightCauchyGreen, deviatoricReferenceStress, pressure ) );

            //Compute the l2norm of the deviatoric stress
            dimVector normDevStress;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress ) );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        }

        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor &stressMeasure,
                                                               const dimVector &cohesion,
                                                               const secondOrderTensor &precedingDeformationGradient,
                                                               const parameterType &frictionAngle, const parameterType &beta,
                                                               dimVector &yieldValue, fourthOrderTensor &dFdStress, secondOrderTensor &dFdc,
                                                               thirdOrderTensor &dFdPrecedingF ){
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
            secondOrderTensor rightCauchyGreen;
            fourthOrderTensor dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            dimVector pressure;
            thirdOrderTensor deviatoricReferenceStress;

            sixthOrderTensor  dDevStressdStress;
            fifthOrderTensor  dDevStressdRCG;
            fourthOrderTensor dPressuredStress;
            thirdOrderTensor  dPressuredRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );

            fifthOrderTensor dDevStressdPrecedingF( tot_dim * sot_dim, 0. );
            thirdOrderTensor dPressuredPrecedingF( dim * sot_dim, 0 );

            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > dRCGdPrecedingF_mat( dRCGdPrecedingF.data( ), sot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > dDevStressdRCG_mat( dDevStressdRCG.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > dDevStressdPrecedingF_mat( dDevStressdPrecedingF.data( ), tot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dPressuredRCG_mat( dPressuredRCG.data( ), dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dPressuredPrecedingF_mat( dPressuredPrecedingF.data( ), dim, sot_dim );

            dDevStressdPrecedingF_mat = ( dDevStressdRCG_mat * dRCGdPrecedingF_mat ).eval( );
            dPressuredPrecedingF_mat = ( dPressuredRCG_mat * dRCGdPrecedingF_mat ).eval( );

            //Compute the l2norm of the deviatoric stress
            dimVector normDevStress;
            fourthOrderTensor dNormDevStressdDevStress;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress ) );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Construct the Jacobians
            dFdStress = BAngle * dPressuredStress;

            Eigen::Map< Eigen::Matrix< floatType, dim, tot_dim, Eigen::RowMajor > > dNormDevStressdDevStress_mat( dNormDevStressdDevStress.data( ), dim, tot_dim );

            Eigen::Map< Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > dDevStressdStress_mat( dDevStressdStress.data( ), tot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< floatType,     dim, tot_dim, Eigen::RowMajor > > dFdStress_mat( dFdStress.data( ), dim, tot_dim );

            dFdStress_mat  = ( dFdStress_mat + dNormDevStressdDevStress_mat * dDevStressdStress_mat ).eval( );

            dFdc = secondOrderTensor( cohesion.size( ) * cohesion.size( ), 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ dFdc[ dim * i + i ] = -AAngle; }

            dFdPrecedingF  = BAngle * dPressuredPrecedingF;

            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dFdPrecedingF_mat( dFdPrecedingF.data( ), dim, sot_dim );

            dFdPrecedingF_mat = ( dFdPrecedingF_mat + dNormDevStressdDevStress_mat * dDevStressdPrecedingF_mat ).eval( );

        }

        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor &stressMeasure,
                                                           const dimVector &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           dimVector &yieldValue, fourthOrderTensor &dFdStress, secondOrderTensor &dFdc,
                                                           thirdOrderTensor &dFdPrecedingF, seventhOrderTensor &d2FdStress2,
                                                           sixthOrderTensor &d2FdStressdPrecedingF ){
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
            constexpr unsigned int fot_dim = tot_dim * dim;
            constexpr unsigned int siot_dim = tot_dim * tot_dim;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            secondOrderTensor rightCauchyGreen;
            fourthOrderTensor dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            dimVector pressure;
            thirdOrderTensor deviatoricReferenceStress;

            sixthOrderTensor dDevStressdStress;
            fifthOrderTensor dDevStressdRCG;
            fourthOrderTensor dPressuredStress;
            thirdOrderTensor dPressuredRCG;

            seventhOrderTensor d2DevStressdStressdRCG;
            fifthOrderTensor d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            fifthOrderTensor dDevStressdPrecedingF( tot_dim * sot_dim, 0. );
            thirdOrderTensor dPressuredPrecedingF( dim * sot_dim, 0 );

            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > dRCGdPrecedingF_mat( dRCGdPrecedingF.data( ), sot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > dDevStressdRCG_mat( dDevStressdRCG.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > dDevStressdPrecedingF_mat( dDevStressdPrecedingF.data( ), tot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dPressuredRCG_mat( dPressuredRCG.data( ), dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dPressuredPrecedingF_mat( dPressuredPrecedingF.data( ), dim, sot_dim );

            dDevStressdPrecedingF_mat = ( dDevStressdRCG_mat * dRCGdPrecedingF_mat ).eval( );
            dPressuredPrecedingF_mat = ( dPressuredRCG_mat * dRCGdPrecedingF_mat ).eval( );

            eighthOrderTensor d2DevStressdStressdPrecedingF( siot_dim * sot_dim, 0 );

            sixthOrderTensor d2PressuredStressdPrecedingF( dim * tot_dim * sot_dim, 0 );

            Eigen::Map< Eigen::Matrix< floatType, siot_dim, sot_dim, Eigen::RowMajor > > d2DevStressdStressdRCG_mat( d2DevStressdStressdRCG.data( ), siot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, siot_dim, sot_dim, Eigen::RowMajor > > d2DevStressdStressdPrecedingF_mat( d2DevStressdStressdPrecedingF.data( ), siot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType,  fot_dim, sot_dim, Eigen::RowMajor > > d2PressuredStressdRCG_mat( d2PressuredStressdRCG.data( ), fot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< floatType,  fot_dim, sot_dim, Eigen::RowMajor > > d2PressuredStressdPrecedingF_mat( d2PressuredStressdPrecedingF.data( ), fot_dim, sot_dim );

            d2DevStressdStressdPrecedingF_mat = ( d2DevStressdStressdRCG_mat * dRCGdPrecedingF_mat ).eval( );
            d2PressuredStressdPrecedingF_mat = ( d2PressuredStressdRCG_mat * dRCGdPrecedingF_mat ).eval( );

            //Compute the l2norm of the deviatoric stress
            thirdOrderTensor normDevStress;
            fourthOrderTensor dNormDevStressdDevStress;
            seventhOrderTensor d2NormDevStressdDevStress2;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
                                                                                                                  dNormDevStressdDevStress,
                                                                                                                  d2NormDevStressdDevStress2 ) )

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Construct the Jacobians
            dFdStress = BAngle * dPressuredStress;

            Eigen::Map< Eigen::Matrix< floatType, dim, tot_dim, Eigen::RowMajor > > dNormDevStressdDevStress_mat( dNormDevStressdDevStress.data( ), dim, tot_dim );

            Eigen::Map< Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > dDevStressdStress_mat( dDevStressdStress.data( ), tot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< floatType,     dim, tot_dim, Eigen::RowMajor > > dFdStress_mat( dFdStress.data( ), dim, tot_dim );

            dFdStress_mat  = ( dFdStress_mat + dNormDevStressdDevStress_mat * dDevStressdStress_mat ).eval( );
    
            dFdc = secondOrderTensor( cohesion.size( ) * cohesion.size( ), 0 );
            for ( unsigned int i = 0; i < dim; i++ ){ dFdc[ dim * i + i ] = -AAngle; }
    
            dFdPrecedingF  = BAngle * dPressuredPrecedingF;

            Eigen::Map< Eigen::Matrix< floatType,     dim, sot_dim, Eigen::RowMajor > > dFdPrecedingF_mat( dFdPrecedingF.data( ), dim, sot_dim );

            dFdPrecedingF_mat = ( dFdPrecedingF_mat + dNormDevStressdDevStress_mat * dDevStressdPrecedingF_mat ).eval( );
    
            //Construct the second-order jacobians
            d2FdStress2 = seventhOrderTensor( dim * tot_dim * tot_dim, 0 );
            d2FdStressdPrecedingF = sixthOrderTensor( dim * tot_dim * sot_dim, 0 );

            sixthOrderTensor   temp_siot1( dim * tot_dim * sot_dim, 0 );
            seventhOrderTensor temp_seot1( dim * tot_dim * tot_dim, 0 );

            for ( unsigned int KLMN = 0; KLMN < fot_dim; KLMN++ ){
                for ( unsigned int ABC = 0; ABC < tot_dim; ABC++ ){
                    for ( unsigned int OP = 0; OP < sot_dim; OP++ ){
                        temp_siot1[ sot_dim * KLMN + OP ]
                            += d2NormDevStressdDevStress2[ tot_dim * KLMN + ABC ]
                             * dDevStressdPrecedingF[ dim * dim * ABC + OP ];
                    }
                }
            }
            for ( unsigned int K = 0; K < 3; K++ ){
                for ( unsigned int ABC = 0; ABC < tot_dim; ABC++ ){
                    for ( unsigned int LMN = 0; LMN < tot_dim; LMN++ ){
                        for ( unsigned int OP = 0; OP < sot_dim; OP++ ){
                            d2FdStressdPrecedingF[ tot_dim * sot_dim * K + sot_dim * LMN + OP ]
                                += temp_siot1[ tot_dim * sot_dim * K + sot_dim * ABC + OP ]
                                 * dDevStressdStress[ tot_dim * ABC + LMN ];
                        }
                    }
                }
            }

            for ( unsigned int K = 0; K < 3; K++ ){
                for ( unsigned int ABC = 0; ABC < tot_dim; ABC++ ){
                    for ( unsigned int LMN = 0; LMN < tot_dim; LMN++ ){
                        for ( unsigned int OPQ = 0; OPQ < tot_dim; OPQ++ ){
                            temp_seot1[ tot_dim * tot_dim * K + tot_dim * LMN + OPQ ]
                                += d2NormDevStressdDevStress2[ tot_dim * tot_dim * K + tot_dim * ABC + OPQ ]
                                 * dDevStressdStress[ tot_dim * ABC + LMN ];
                        }
                    }
                }
            }

            for ( unsigned int K = 0; K < 3; K++ ){
                for ( unsigned int LMN = 0; LMN < tot_dim; LMN++ ){
                    for ( unsigned int ABC = 0; ABC < tot_dim; ABC++ ){
                        for ( unsigned int OPQ = 0; OPQ < tot_dim; OPQ++ ){
                            d2FdStress2[ tot_dim * tot_dim * K + tot_dim * LMN + OPQ ]
                                += temp_seot1[ tot_dim * tot_dim * K + tot_dim * LMN + ABC ]
                                 * dDevStressdStress[ tot_dim * ABC + OPQ ];
                        }
                    }
                }
            }

            d2FdStressdPrecedingF += BAngle * d2PressuredStressdPrecedingF;

            Eigen::Map< Eigen::Matrix< floatType,     dim, tot_dim * sot_dim, Eigen::RowMajor > > d2FdStressdPrecedingF_mat( d2FdStressdPrecedingF.data( ), dim, tot_dim * sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, tot_dim * sot_dim, Eigen::RowMajor > > d2DevStressdStressdPrecedingF_mat2( d2DevStressdStressdPrecedingF.data( ), tot_dim, tot_dim * sot_dim );

            d2FdStressdPrecedingF_mat = ( d2FdStressdPrecedingF_mat + dNormDevStressdDevStress_mat * d2DevStressdStressdPrecedingF_mat2 ).eval( );

        }

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient ){
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

            TARDIGRADE_ERROR_TOOLS_CHECK( inverseElasticRightCauchyGreen.size() == dim * dim, "The inverse elastic right Cauchy-Green deformation tensor must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( macroFlowDirection.size() == dim * dim, "The macro flow direction tensor must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( microFlowDirection.size() == dim * dim, "The micro flow direction tensor must be 3D" );

            //Compute the macro-scale velocity gradient
            plasticMacroVelocityGradient = secondOrderTensor( dim * dim, 0 );

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
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient,
                                                  secondOrderTensor &dPlasticMacroLdMacroGamma,
                                                  secondOrderTensor &dPlasticMacroLdMicroGamma ){
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

            dPlasticMacroLdMacroGamma = secondOrderTensor( dim * dim, 0 );
            dPlasticMacroLdMicroGamma = secondOrderTensor( dim * dim, 0 );

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
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient,
                                                  secondOrderTensor &dPlasticMacroLdMacroGamma,
                                                  secondOrderTensor &dPlasticMacroLdMicroGamma,
                                                  fourthOrderTensor &dPlasticMacroLdElasticRCG,
                                                  fourthOrderTensor &dPlasticMacroLdMacroFlowDirection,
                                                  fourthOrderTensor &dPlasticMacroLdMicroFlowDirection ){
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

            dPlasticMacroLdElasticRCG = fourthOrderTensor( sot_dim * sot_dim, 0 );
            dPlasticMacroLdMacroFlowDirection = fourthOrderTensor( sot_dim * sot_dim, 0 );
            dPlasticMacroLdMicroFlowDirection = fourthOrderTensor( sot_dim * sot_dim, 0 );

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

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient ){
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

            plasticMicroVelocityGradient = secondOrderTensor( sot_dim, 0 );

            //NOTE: I'm making the second inverse elastic Psi be the transpose of what was done previously.
            //      I think the way it was is a bug since it isn't consistent with the form in my dissertation.
            secondOrderTensor temp_sot( sot_dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        plasticMicroVelocityGradient[ dim * Bb + Kb ] += microGamma * inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Kb + Lb ];

                    }

                }

            }

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        temp_sot[ dim * Bb + Kb ]
                            += plasticMicroVelocityGradient[ dim * Bb + Lb ] * inverseElasticPsi[ dim * Kb + Lb ];

                    }

                }

            }

            std::fill( plasticMicroVelocityGradient.begin( ),
                       plasticMicroVelocityGradient.end( ),
                       0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        plasticMicroVelocityGradient[ dim * Bb + Kb ]
                            += temp_sot[ dim * Bb + Lb ] * elasticMicroRightCauchyGreen[ dim * Lb + Kb ];

                    }

                }

            }

        }

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient,
                                                  secondOrderTensor &dPlasticMicroLdMicroGamma ){
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
            constexpr unsigned int sot_dim = dim * dim;

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticMicroRightCauchyGreen.size() == dim * dim, "The elastic micro right Cauchy-Green deformation tensor is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticPsi.size() == dim * dim, "The elastic micro deformation tensor Psi is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( inverseElasticPsi.size() == dim * dim, "The inverse of the elastic micro deformation tensor Psi is not 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( microFlowDirection.size() == dim * dim, "The micro flow direction of the elastic micro plastic flow direction is not 3D" );

            plasticMicroVelocityGradient = secondOrderTensor( dim * dim, 0 );
            dPlasticMicroLdMicroGamma = secondOrderTensor( dim * dim, 0 );

            secondOrderTensor temp_sot( sot_dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        dPlasticMicroLdMicroGamma[ dim * Bb + Kb ] += inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Kb + Lb ];

                    }

                }

            }

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        temp_sot[ dim * Bb + Kb ]
                            += dPlasticMicroLdMicroGamma[ dim * Bb + Lb ] * inverseElasticPsi[ dim * Kb + Lb ];

                    }

                }

            }

            std::fill( dPlasticMicroLdMicroGamma.begin( ),
                       dPlasticMicroLdMicroGamma.end( ),
                       0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        dPlasticMicroLdMicroGamma[ dim * Bb + Kb ]
                            += temp_sot[ dim * Bb + Lb ] * elasticMicroRightCauchyGreen[ dim * Lb + Kb ];

                    }

                }

            }

            plasticMicroVelocityGradient = microGamma * dPlasticMicroLdMicroGamma;

        }

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient,
                                                  secondOrderTensor &dPlasticMicroLdMicroGamma,
                                                  fourthOrderTensor &dPlasticMicroLdElasticMicroRCG,
                                                  fourthOrderTensor &dPlasticMicroLdElasticPsi,
                                                  fourthOrderTensor &dPlasticMicroLdMicroFlowDirection ){
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
            constexpr unsigned int fot_dim = sot_dim * sot_dim;

            TARDIGRADE_ERROR_TOOLS_CATCH(

                computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen,
                                                     elasticPsi, inverseElasticPsi, microFlowDirection,
                                                     plasticMicroVelocityGradient, dPlasticMicroLdMicroGamma );

            )

            //Assemble the Jacobians
            dPlasticMicroLdElasticMicroRCG = fourthOrderTensor( sot_dim * sot_dim, 0 );

            dPlasticMicroLdElasticPsi = fourthOrderTensor( sot_dim * sot_dim, 0 );

            dPlasticMicroLdMicroFlowDirection = fourthOrderTensor( sot_dim * sot_dim, 0 );

            secondOrderTensor temp_sot1( sot_dim, 0 );

            secondOrderTensor temp_sot1a( sot_dim, 0 );

            secondOrderTensor temp_sot1b( sot_dim, 0 );

            for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        temp_sot1[ dim * Bb + Kb ]  -= microGamma * inverseElasticPsi[ dim * Ob + Kb ] * elasticMicroRightCauchyGreen[ dim * Ob + Bb ];

                        temp_sot1a[ dim * Bb + Kb ] += microFlowDirection[ dim * Ob + Bb ] * inverseElasticPsi[ dim * Kb + Ob ];

                    }

                }

            }

            for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        temp_sot1b[ dim * Bb + Kb ] += inverseElasticPsi[ dim * Bb + Ob ] * temp_sot1a[ dim * Ob + Kb ];

                    }

                }

            }

            fourthOrderTensor temp_fot1( fot_dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                        dPlasticMicroLdElasticMicroRCG[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Kb ]
                            += microGamma * temp_sot1b[ dim * Bb + Ob ];

                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMicroLdElasticPsi[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
                                += temp_sot1[ dim * Kb + Ob ] * temp_sot1b[ dim * Bb + Pb ]
                                 - inverseElasticPsi[ dim * Bb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Kb ];

                            dPlasticMicroLdMicroFlowDirection[ dim * sot_dim * Bb + sot_dim * Kb + dim * Ob + Pb ]
                                -= inverseElasticPsi[ dim * Bb + Pb ] * temp_sot1[ dim * Kb + Ob ];

                        }

                    }

                }

            }

        }

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient ){
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

            thirdOrderTensor skewTerm;
            return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                elasticGamma, microGradientFlowDirection,
                                                                plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                                skewTerm );

        }

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          thirdOrderTensor &skewTerm ){
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
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;

            TARDIGRADE_ERROR_TOOLS_CHECK( microGradientGamma.size() == dim, "The micro gradient plastic multiplier must have a length of 3" );

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticPsi.size() == sot_dim, "The elastic micro deformation measure Psi must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( inverseElasticPsi.size() == sot_dim, "The inverse elastic micro deformation measure Psi must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( elasticGamma.size() == tot_dim, "The elastic higher order deformation measure Gamma must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( microGradientFlowDirection.size() == dim * tot_dim, "The micro gradient flow direction must be 3D" );

            TARDIGRADE_ERROR_TOOLS_CHECK( plasticMicroVelocityGradient.size() == sot_dim, "The plastic micro velocity gradient must be 3D" );

            //Assemble the 'skew' term
            skewTerm = thirdOrderTensor( tot_dim, 0 );
            secondOrderTensor temp_sot( sot_dim, 0 );
            thirdOrderTensor temp_tot( tot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        temp_sot[ dim * Db + Kb ] += plasticMicroVelocityGradient[ dim * Db + Mb ]
                                                   * inverseElasticPsi[ dim * Mb + Kb ];

                        for ( unsigned int Cb = 0; Cb < dim; Cb++ ){

                            temp_tot[ dim * dim * Db + dim * Kb + Cb ] -= inverseElasticPsi[ dim * Db + Mb ] * elasticGamma[ dim * dim * Mb + dim * Kb + Cb ];

                        }

                    }

                }

            }

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Cb = 0; Cb < dim; Cb++ ){

                            skewTerm[ dim * dim * Db + dim * Mb + Kb ]
                                += temp_sot[ dim * Db + Cb ] * elasticGamma[ dim * dim * Cb + dim * Mb + Kb ]
                                 + plasticMicroVelocityGradient[ dim * Cb + Mb ] * temp_tot[ dim * dim * Db + dim * Cb + Kb ];

                        }

                    }

                }

            }

            plasticMicroGradientVelocityGradient = thirdOrderTensor( dim * dim * dim, 0 );

            std::fill( temp_tot.begin( ), temp_tot.end( ), 0 );

            for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                            temp_tot[ dim * dim * Nb + dim * Mb + Kb ]
                                += microGradientGamma[ Lb ] * microGradientFlowDirection[ dim * dim * dim * Lb + dim * dim * Mb + dim * Kb + Nb ]
                                 + elasticPsi[ dim * Kb + Lb ] * skewTerm[ dim * dim * Lb + dim * Nb + Mb ];

                        }

                    }

                }

            }

            for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                            plasticMicroGradientVelocityGradient[ dim * dim * Nb + dim * Mb + Kb ]
                                += inverseElasticPsi[ dim * Nb + Lb ] * temp_tot[ dim * dim * Mb + dim * Kb + Lb ];

                        }

                    }

                }

            }

        }

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL ){
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

            thirdOrderTensor skewTerm;
            return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi, elasticGamma,
                                                                microGradientFlowDirection, plasticMicroVelocityGradient,
                                                                plasticMicroGradientVelocityGradient, skewTerm,
                                                                dPlasticMicroGradientLdMicroGradientGamma,
                                                                dPlasticMicroGradientLdPlasticMicroL );
        }

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          thirdOrderTensor &skewTerm,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL ){
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

            dPlasticMicroGradientLdPlasticMicroL = fifthOrderTensor( tot_dim * sot_dim, 0 );
            dPlasticMicroGradientLdMicroGradientGamma = fourthOrderTensor( tot_dim * dim, 0 );

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

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor  &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor  &dPlasticMicroGradientLdPlasticMicroL,
                                                          fifthOrderTensor  &dPlasticMicroGradientLdElasticPsi,
                                                          sixthOrderTensor  &dPlasticMicroGradientLdElasticGamma,
                                                          sixthOrderTensor  &dPlasticMicroGradientLdMicroGradientFlowDirection ){
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

            thirdOrderTensor skewTerm;
            TARDIGRADE_ERROR_TOOLS_CATCH(

                computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection,
                                                             plasticMicroVelocityGradient,
                                                             plasticMicroGradientVelocityGradient, skewTerm,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL );

            )

            dPlasticMicroGradientLdElasticPsi   = fifthOrderTensor( tot_dim * sot_dim, 0 );

            dPlasticMicroGradientLdElasticGamma = sixthOrderTensor( tot_dim * tot_dim, 0 );

            dPlasticMicroGradientLdMicroGradientFlowDirection = seventhOrderTensor( tot_dim * dim * tot_dim, 0 );

            secondOrderTensor temp_sot1( sot_dim, 0 );

            thirdOrderTensor temp_tot1( tot_dim, 0 );

            thirdOrderTensor temp_tot2( tot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        temp_sot1[ dim * Db + Kb ] += plasticMicroVelocityGradient[ dim * Db + Mb ] * inverseElasticPsi[ dim * Mb + Kb ];

                        for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                            temp_tot1[ dim * dim * Db + dim * Kb + Ob ] += inverseElasticPsi[ dim * Db + Mb ]
                                                                         * elasticGamma[ dim * dim * Mb + dim * Kb + Ob ]; 

                        }

                    }

                }

            }

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                            temp_tot2[ dim * dim * Db + dim * Kb + Ob ] += plasticMicroVelocityGradient[ dim * Mb + Kb ] * temp_tot1[ dim * dim * Db + dim * Mb + Ob ];



                        }

                    }

                }

            }

            for ( unsigned int Db = 0; Db < dim; Db++ ){

                for ( unsigned int Mb = 0; Mb < dim; Mb++ ){

                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                        for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                            dPlasticMicroGradientLdElasticGamma[ dim * dim * tot_dim * Db + dim * tot_dim * Mb + tot_dim * Kb + dim * dim * Ob + dim * Mb + Kb ]
                                += temp_sot1[ dim * Db + Ob ];

                            for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                                dPlasticMicroGradientLdElasticPsi[ dim * dim * sot_dim * Db + dim * sot_dim * Mb + sot_dim * Kb + dim * Ob + Pb ]
                                    += inverseElasticPsi[ dim * Db + Ob ]
                                     * ( skewTerm[ dim * dim * Pb + dim * Mb + Kb ] - plasticMicroGradientVelocityGradient[ dim * dim * Pb + dim * Mb + Kb ] );

                                dPlasticMicroGradientLdElasticGamma[ dim * dim * tot_dim * Db + dim * tot_dim * Mb + tot_dim * Kb + dim * dim * Ob + dim * Pb + Kb ]
                                    -= plasticMicroVelocityGradient[ dim * Pb + Mb ] * inverseElasticPsi[ dim * Db + Ob ];

                                dPlasticMicroGradientLdMicroGradientFlowDirection[ dim * dim * tot_dim * dim * Db + dim * tot_dim * dim * Mb + tot_dim * dim * Kb + dim * dim * dim * Ob + dim * dim * Kb + dim * Pb + Mb ]
                                    += inverseElasticPsi[ dim * Db + Pb ] * microGradientGamma[ Ob ];

                                dPlasticMicroGradientLdElasticPsi[ dim * dim * sot_dim * Db + dim * sot_dim * Mb + sot_dim * Kb + dim * Ob + Pb ]
                                    += inverseElasticPsi[ dim * Db + Ob ] * temp_tot2[ dim * dim * Pb + dim * Mb + Kb ]
                                     - temp_sot1[ dim * Db + Ob ] * temp_tot1[ dim * dim * Pb + dim * Mb + Kb ];

                            }

                        }

                    }

                }

            }

        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor  &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor        &currentPlasticMicroGradient,
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

            sixthOrderTensor LHS;

            evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                       currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                       previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                       previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                       previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient, LHS,
                                       alpha );
        }

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const secondOrderTensor &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        sixthOrderTensor &LHS,
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
            thirdOrderTensor DtAtilde( tot_dim, 0 );
            fourthOrderTensor previousFourthA( fot_dim, 0 );
            fourthOrderTensor currentFourthA( fot_dim, 0 );

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
            thirdOrderTensor RHS = DtAtilde;
            RHS += previousPlasticMicroGradient;
            LHS = sixthOrderTensor( tot_dim * tot_dim, 0 );
            for ( unsigned int i = 0; i < tot_dim; i++ ){ LHS[ tot_dim * i + i ] = 1; }

            for ( unsigned int Db = 0; Db < dim; Db++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
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
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor  &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor        &currentPlasticMicroGradient,
                                        fifthOrderTensor        &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        fifthOrderTensor        &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        fifthOrderTensor        &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        sixthOrderTensor        &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
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
            sixthOrderTensor LHS;
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
            sixthOrderTensor   negdRdCurrentDtAtilde( dim * dim * dim * dim * dim * dim, 0 );
            seventhOrderTensor negdRdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim, 0 );

            //Also assemble jacobians of the A terms
            fifthOrderTensor dCurrentDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            sixthOrderTensor dCurrentDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );
            sixthOrderTensor dCurrentFourthAdMacroVelocityGradient( fot_dim * sot_dim, 0 );
            sixthOrderTensor dCurrentFourthAdMicroVelocityGradient( fot_dim * sot_dim, 0 );

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
            sixthOrderTensor dCurrentPlasticMicroGradientdCurrentDTAtilde( dim * dim * dim * dim * dim * dim );
            seventhOrderTensor dCurrentPlasticMicroGradientdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim );

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
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor  &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient,
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

            //Compute the new currentPlasticMicroGradient
            sixthOrderTensor LHS;
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
            sixthOrderTensor   negdRdCurrentDtAtilde( tot_dim * tot_dim, 0 );
            seventhOrderTensor negdRdCurrentFourthA( tot_dim * fot_dim, 0 );

            //Also assemble jacobians of the A terms
            fifthOrderTensor dCurrentDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            sixthOrderTensor dCurrentDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );
            sixthOrderTensor dCurrentFourthAdMacroVelocityGradient( fot_dim * sot_dim, 0 );
            sixthOrderTensor dCurrentFourthAdMicroVelocityGradient( fot_dim * sot_dim, 0 );

            fifthOrderTensor dPreviousDTAtildedPlasticMicroDeformation( tot_dim * sot_dim, 0 );
            sixthOrderTensor dPreviousDTAtildedPlasticMicroGradientVelocityGradient( tot_dim * tot_dim, 0 );

            sixthOrderTensor dRHSdPreviousPlasticMicroGradient( tot_dim * tot_dim, 0 );
            fifthOrderTensor dRHSdPreviousPlasticMacroVelocityGradient( tot_dim * sot_dim, 0 );
            fifthOrderTensor dRHSdPreviousPlasticMicroVelocityGradient( tot_dim * sot_dim, 0 );

            for ( unsigned int Db = 0; Db < dim; Db++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                        negdRdCurrentDtAtilde[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Db + dim * B + Kb ] += 1;

                        dCurrentFourthAdMicroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * B + dim * sot_dim * Kb + sot_dim * Kb + dim * Db + B ] += 1;

                        dCurrentFourthAdMacroVelocityGradient[ dim * dim * dim * sot_dim * Db + dim * dim * sot_dim * Db + dim * sot_dim * B + sot_dim * Kb + dim * Kb + B ] -= 1;

                        dRHSdPreviousPlasticMicroGradient[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Db + dim * B + Kb ] += 1;

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

                            dRHSdPreviousPlasticMicroGradient[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Rb + dim * B + Kb ]
                                += Dt * ( 1 - alpha ) * previousPlasticMicroVelocityGradient[ dim * Db + Rb ];

                            dRHSdPreviousPlasticMicroGradient[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Db + dim * B + Rb ]
                                -= Dt * ( 1 - alpha ) * previousPlasticMacroVelocityGradient[ dim * Rb + Kb ];

                            for ( unsigned int S = 0; S < dim; S++ ){

                                negdRdCurrentFourthA[ dim * dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * dim * B + dim * dim * dim * dim * Kb + dim * dim * dim * Db + dim * dim * Rb + dim * Kb + S ]
                                    += Dt * alpha * currentPlasticMicroGradient[ dim * dim * Rb + dim * B + S ];

                            }
                        }
                    }
                }
            }

            sixthOrderTensor dCurrentPlasticMicroGradientdCurrentDTAtilde( tot_dim * tot_dim, 0 );
            seventhOrderTensor dCurrentPlasticMicroGradientdCurrentFourthA( tot_dim * fot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient = sixthOrderTensor( tot_dim * tot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient = fifthOrderTensor( tot_dim * sot_dim, 0 );
            dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient = fifthOrderTensor( tot_dim * sot_dim, 0 );

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
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor  &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor  &currentPlasticMicroGradient,
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
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
//                                                                    currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
//                                                                    alphaMacro );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            1. - alphaMicro, 1 );
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
//                                                                    currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
//                                                                    alphaMicro );
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
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor  &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor  &currentPlasticMicroGradient,
                                       fourthOrderTensor &dPlasticFdPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPlasticMicroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPlasticMacroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPlasticMicroL,
                                       sixthOrderTensor  &dPlasticMicroGradientdPlasticMicroGradientL,
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
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
//                                                                    currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
//                                                                    dPlasticFdPlasticMacroL, alphaMacro );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, 1. - alphaMicro, 1 );
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
//                                                                    currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
//                                                                    dPlasticMicroDeformationdPlasticMicroL, alphaMicro );
            )

            fifthOrderTensor dPlasticMicroGradientdPlasticMicroDeformation;
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
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor  &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor  &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor  &currentPlasticMicroGradient,
                                       fourthOrderTensor &dPlasticFdPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPlasticMicroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPlasticMacroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPlasticMicroL,
                                       sixthOrderTensor  &dPlasticMicroGradientdPlasticMicroGradientL,
                                       fourthOrderTensor &dPlasticFdPreviousPlasticF,
                                       fourthOrderTensor &dPlasticFdPreviousPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                       fourthOrderTensor &dPlasticMicroDeformationdPreviousPlasticMicroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                       fifthOrderTensor  &dPlasticMicroGradientdPreviousPlasticMicroGradient,
                                       fifthOrderTensor  &dPlasticMicroGradientdPreviousPlasticMacroL,
                                       fifthOrderTensor  &dPlasticMicroGradientdPreviousPlasticMicroL,
                                       sixthOrderTensor  &dPlasticMicroGradientdPreviousPlasticMicroGradientL,
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
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
//                                                                    currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
//                                                                    dPlasticFdPlasticMacroL, dPlasticFdPreviousPlasticF, dPlasticFdPreviousPlasticMacroL, alphaMacro );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeConstitutiveTools::evolveFFlatJ( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPreviousPlasticMicroL, 1. - alphaMicro, 1 );
//                tardigradeConstitutiveTools::evolveFExponentialMap( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
//                                                                    currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
//                                                                    dPlasticMicroDeformationdPlasticMicroL, dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
//                                                                    dPlasticMicroDeformationdPreviousPlasticMicroL, alphaMicro );
            )

            fifthOrderTensor dPlasticMicroGradientdPlasticMicroDeformation;
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

            secondOrderTensor Fp;

            secondOrderTensor chip;

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
            secondOrderTensor PK2Stress(                     stress->begin( ),               stress->begin( ) + 1 * sot_dim );

            secondOrderTensor referenceSymmetricMicroStress( stress->begin( ) + 1 * sot_dim, stress->begin( ) + 2 * sot_dim );

            thirdOrderTensor referenceHigherOrderStress(    stress->begin( ) + 2 * sot_dim, stress->begin( ) + 2 * sot_dim + tot_dim );

            // Push the stresses forward to the current configuration of the plastic configuration
            secondOrderTensor macroDrivingStress;

            secondOrderTensor symmetricMicroDrivingStress;

            thirdOrderTensor higherOrderDrivingStress;

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

            const fourthOrderTensor *dF1dF;

            const floatVector *dF1dFn;

            floatVector dChipdSubChis;

            const fourthOrderTensor *dChi1dChi;

            const floatVector *dChi1dChin;

            secondOrderTensor Fp;

            secondOrderTensor chip;

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
            fourthOrderTensor dFpdF(  sot_dim * sot_dim, 0 );

            floatVector dFpdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            fourthOrderTensor dChipdChi(  sot_dim * sot_dim, 0 );

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
            secondOrderTensor PK2Stress(                     stress->begin( ),               stress->begin( ) + 1 * sot_dim );

            secondOrderTensor referenceSymmetricMicroStress( stress->begin( ) + 1 * sot_dim, stress->begin( ) + 2 * sot_dim );;

            thirdOrderTensor  referenceHigherOrderStress(    stress->begin( ) + 2 * sot_dim, stress->begin( ) + 2 * sot_dim + tot_dim );

            // Push the stresses forward to the current configuration of the plastic configuration
            secondOrderTensor macroDrivingStress;

            secondOrderTensor symmetricMicroDrivingStress;

            thirdOrderTensor higherOrderDrivingStress;

            fourthOrderTensor dMacrodFp;

            fourthOrderTensor dMacrodPK2;

            fourthOrderTensor dMicrodFp;

            fourthOrderTensor dMicrodSigma;

            fifthOrderTensor  dHigherdFp;

            fifthOrderTensor  dHigherdChip;

            sixthOrderTensor  dHigherdM;

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

            const dimVector  *microGradientCohesion;

            const secondOrderTensor *macroDrivingStress;

            const secondOrderTensor *microDrivingStress;

            const thirdOrderTensor  *microGradientDrivingStress;

            const secondOrderTensor *precedingDeformationGradient;

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

            secondOrderTensor dMacroFlowdDrivingStress, dMicroFlowdDrivingStress,
                              dMacroFlowdPrecedingF,    dMicroFlowdPrecedingF;

            fourthOrderTensor dMicroGradientFlowdDrivingStress;
            secondOrderTensor dMicroGradientFlowdCohesion;
            thirdOrderTensor  dMicroGradientFlowdPrecedingF;

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

            const dimVector *microGradientCohesion;

            const secondOrderTensor *macroDrivingStress;

            const secondOrderTensor *microDrivingStress;

            const thirdOrderTensor  *microGradientDrivingStress;

            const fourthOrderTensor *dMacroDrivingStressdStress;

            const fourthOrderTensor *dMacroDrivingStressdF;

            const floatVector *dMacroDrivingStressdFn;

            const fourthOrderTensor *dMicroDrivingStressdStress;

            const fourthOrderTensor *dMicroDrivingStressdF;

            const floatVector *dMicroDrivingStressdFn;

            const sixthOrderTensor *dMicroGradientDrivingStressdStress;

            const fifthOrderTensor *dMicroGradientDrivingStressdF;

            const floatVector *dMicroGradientDrivingStressdFn;

            const fifthOrderTensor *dMicroGradientDrivingStressdChi;

            const floatVector *dMicroGradientDrivingStressdChin;

            const secondOrderTensor *precedingDeformationGradient;

            const fourthOrderTensor *dPrecedingFdF;

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

            secondOrderTensor dMacroFlowdDrivingStress, dMicroFlowdDrivingStress,
                              dMacroFlowdPrecedingF,    dMicroFlowdPrecedingF;

            fourthOrderTensor dMicroGradientFlowdDrivingStress;
            secondOrderTensor dMicroGradientFlowdCohesion;
            thirdOrderTensor  dMicroGradientFlowdPrecedingF;

            fourthOrderTensor d2MacroFlowdDrivingStress2,         d2MacroFlowdDrivingStressdPrecedingF,
                              d2MicroFlowdDrivingStress2,         d2MicroFlowdDrivingStressdPrecedingF;

            seventhOrderTensor d2MicroGradientFlowdDrivingStress2;
            sixthOrderTensor   d2MicroGradientFlowdDrivingStressdPrecedingF;

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

            fourthOrderTensor  d2MacroFlowdDrivingStressdMacroStress( sot_dim * sot_dim, 0 );

            fourthOrderTensor  d2MicroFlowdDrivingStressdMicroStress( sot_dim * sot_dim, 0 );

            seventhOrderTensor d2MicroGradientFlowdDrivingStressdMicroGradientStress( dim * tot_dim * tot_dim, 0 );

            fourthOrderTensor  d2MacroFlowdDrivingStressdF( sot_dim * sot_dim, 0 );

            fourthOrderTensor  d2MicroFlowdDrivingStressdF( sot_dim * sot_dim, 0 );

            sixthOrderTensor   d2MicroGradientFlowdDrivingStressdF( dim * tot_dim * sot_dim, 0 );

            sixthOrderTensor   d2MicroGradientFlowdDrivingStressdChi( dim * tot_dim * sot_dim, 0 );

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

            TARDIGRADE_ERROR_TOOLS_CHECK( get_macroHardeningParameters( )->size( ) == 2, "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_macroHardeningParameters( )->size( ) ) );

            TARDIGRADE_ERROR_TOOLS_CHECK( get_microHardeningParameters( )->size( ) == 2, "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microHardeningParameters( )->size( ) ) );

            TARDIGRADE_ERROR_TOOLS_CHECK( get_microGradientHardeningParameters( )->size( ) == 2, "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microGradientHardeningParameters( )->size( ) ) );

            floatType macroCohesion           = ( *get_macroHardeningParameters( ) )[ 0 ] + ( *get_macroHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 0 ];

            floatType microCohesion           = ( *get_microHardeningParameters( ) )[ 0 ] + ( *get_microHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 1 ];

            dimVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * dimVector( plasticStrainLikeISVs->begin( ) + 2,
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

            dimVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * dimVector( plasticStrainLikeISVs->begin( ) + 2,
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

            const secondOrderTensor *dMicroGradientFlowdc;

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

            const secondOrderTensor *dMicroGradientFlowdc;

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

            const secondOrderTensor *macroDrivingStress;

            const secondOrderTensor *microDrivingStress;

            const thirdOrderTensor  *microGradientDrivingStress;

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const dimVector   *microGradientCohesion;

            const secondOrderTensor *precedingDeformationGradient;

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

            dimVector microGradientYield;

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

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const secondOrderTensor *macroDrivingStress;

            const secondOrderTensor *microDrivingStress;

            const thirdOrderTensor  *microGradientDrivingStress;

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const dimVector   *microGradientCohesion;

            const floatVector *dMacroCohesiondStateVariables;

            const fourthOrderTensor *dMacroDrivingStressdStress;

            const fourthOrderTensor *dMacroDrivingStressdF;

            const floatVector *dMacroDrivingStressdFn;

            const floatVector *dMicroCohesiondStateVariables;

            const fourthOrderTensor *dMicroDrivingStressdStress;

            const fourthOrderTensor *dMicroDrivingStressdF;

            const floatVector *dMicroDrivingStressdFn;

            const floatVector *dMicroGradientCohesiondStateVariables;

            const floatVector *dMicroGradientDrivingStressdStress;

            const fourthOrderTensor *dMicroGradientDrivingStressdF;

            const floatVector *dMicroGradientDrivingStressdFn;

            const fourthOrderTensor *dMicroGradientDrivingStressdChi;

            const floatVector *dMicroGradientDrivingStressdChin;

            const floatVector *precedingDeformationGradient;

            const fourthOrderTensor *dPrecedingFdF;

            const floatVector *dPrecedingFdFn;

            const floatVector *macroYieldParameters = get_macroYieldParameters( );

            const floatVector *microYieldParameters = get_microYieldParameters( );

            const floatVector *microGradientYieldParameters = get_microGradientYieldParameters( );

            setDataStorageBase< floatType > macroYield;

            setDataStorageBase< floatType > microYield;

            setDataStorageBase< dimVector > microGradientYield;

            setDataStorageBase< secondOrderTensor > dMacroYielddStress;

            setDataStorageBase< floatVector > dMacroYielddStateVariables;

            setDataStorageBase< secondOrderTensor > dMacroYielddF;

            setDataStorageBase< floatVector > dMacroYielddFn;

            setDataStorageBase< secondOrderTensor > dMicroYielddStress;

            setDataStorageBase< floatVector > dMicroYielddStateVariables;

            setDataStorageBase< secondOrderTensor > dMicroYielddF;

            setDataStorageBase< floatVector > dMicroYielddFn;

            setDataStorageBase< fourthOrderTensor > dMicroGradientYielddStress;

            setDataStorageBase< floatVector > dMicroGradientYielddStateVariables;

            setDataStorageBase< fourthOrderTensor > dMicroGradientYielddF;

            setDataStorageBase< floatVector > dMicroGradientYielddFn;

            setDataStorageBase< fourthOrderTensor > dMicroGradientYielddChi;

            setDataStorageBase< floatVector > dMicroGradientYielddChin;

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

                macroDrivingStress                    = get_previousMacroDrivingStress( );

                microDrivingStress                    = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress            = get_previousHigherOrderDrivingStress( );

                macroCohesion                         = get_previousMacroCohesion( );

                microCohesion                         = get_previousMicroCohesion( );

                microGradientCohesion                 = get_previousMicroGradientCohesion( );

                macroYield                            = get_setDataStorage_previousMacroYield( );

                microYield                            = get_setDataStorage_previousMicroYield( );

                microGradientYield                    = get_setDataStorage_previousMicroGradientYield( );

                dMacroYielddStress                    = get_setDataStorage_previousdMacroYielddStress( );

                dMacroYielddStateVariables            = get_setDataStorage_previousdMacroYielddStateVariables( );

                dMacroYielddF                         = get_setDataStorage_previousdMacroYielddF( );

                dMacroYielddFn                        = get_setDataStorage_previousdMacroYielddFn( );

                dMicroYielddStress                    = get_setDataStorage_previousdMicroYielddStress( );

                dMicroYielddStateVariables            = get_setDataStorage_previousdMicroYielddStateVariables( );

                dMicroYielddF                         = get_setDataStorage_previousdMicroYielddF( );

                dMicroYielddFn                        = get_setDataStorage_previousdMicroYielddFn( );

                dMicroGradientYielddStress            = get_setDataStorage_previousdMicroGradientYielddStress( );

                dMicroGradientYielddStateVariables    = get_setDataStorage_previousdMicroGradientYielddStateVariables( );

                dMicroGradientYielddF                 = get_setDataStorage_previousdMicroGradientYielddF( );

                dMicroGradientYielddFn                = get_setDataStorage_previousdMicroGradientYielddFn( );

                dMicroGradientYielddChi               = get_setDataStorage_previousdMicroGradientYielddChi( );

                dMicroGradientYielddChin              = get_setDataStorage_previousdMicroGradientYielddChin( );

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

                macroDrivingStress                 = get_macroDrivingStress( );

                microDrivingStress                 = get_symmetricMicroDrivingStress( );

                microGradientDrivingStress         = get_higherOrderDrivingStress( );

                macroCohesion                      = get_macroCohesion( );

                microCohesion                      = get_microCohesion( );

                microGradientCohesion              = get_microGradientCohesion( );

                macroYield                         = get_setDataStorage_macroYield( );

                microYield                         = get_setDataStorage_microYield( );

                microGradientYield                 = get_setDataStorage_microGradientYield( );

                dMacroYielddStress                 = get_setDataStorage_dMacroYielddStress( );

                dMacroYielddStateVariables         = get_setDataStorage_dMacroYielddStateVariables( );

                dMacroYielddF                      = get_setDataStorage_dMacroYielddF( );

                dMacroYielddFn                     = get_setDataStorage_dMacroYielddFn( );

                dMicroYielddStress                 = get_setDataStorage_dMicroYielddStress( );

                dMicroYielddStateVariables         = get_setDataStorage_dMicroYielddStateVariables( );

                dMicroYielddF                      = get_setDataStorage_dMicroYielddF( );

                dMicroYielddFn                     = get_setDataStorage_dMicroYielddFn( );

                dMicroGradientYielddStress         = get_setDataStorage_dMicroGradientYielddStress( );

                dMicroGradientYielddStateVariables = get_setDataStorage_dMicroGradientYielddStateVariables( );

                dMicroGradientYielddF              = get_setDataStorage_dMicroGradientYielddF( );

                dMicroGradientYielddFn             = get_setDataStorage_dMicroGradientYielddFn( );

                dMicroGradientYielddChi            = get_setDataStorage_dMicroGradientYielddChi( );

                dMicroGradientYielddChin           = get_setDataStorage_dMicroGradientYielddChin( );

            }

            secondOrderTensor dMacroYielddDrivingStress;

            floatType   dMacroYielddCohesion;

            secondOrderTensor dMacroYielddPrecedingF;

            secondOrderTensor dMicroYielddDrivingStress;

            floatType   dMicroYielddCohesion;

            secondOrderTensor dMicroYielddPrecedingF;

            fourthOrderTensor dMicroGradientYielddDrivingStress;

            secondOrderTensor dMicroGradientYielddCohesion;

            thirdOrderTensor dMicroGradientYielddPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                                                                                        ( *macroYieldParameters )[ 0 ], ( *macroYieldParameters )[ 1 ],
                                                                                        *macroYield.value, dMacroYielddDrivingStress, dMacroYielddCohesion, dMacroYielddPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                                                                                        ( *microYieldParameters )[ 0 ], ( *microYieldParameters )[ 1 ],
                                                                                        *microYield.value, dMicroYielddDrivingStress, dMicroYielddCohesion, dMicroYielddPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
                                                                                        ( *microGradientYieldParameters )[ 0 ], ( *microGradientYieldParameters )[ 1 ],
                                                                                        *microGradientYield.value, dMicroGradientYielddDrivingStress, dMicroGradientYielddCohesion,
                                                                                        dMicroGradientYielddPrecedingF ) );

            auto map_dMacroYielddDrivingStress         = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dMacroYielddDrivingStress.data( ) );
            auto map_dMacroYielddPrecedingF            = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dMacroYielddPrecedingF.data( ) );
            auto map_dMicroYielddDrivingStress         = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dMicroYielddDrivingStress.data( ) );
            auto map_dMicroYielddPrecedingF            = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dMicroYielddPrecedingF.data( ) );
            auto map_dMicroGradientYielddDrivingStress = getFixedSizeMatrixMap< floatType,       3, tot_dim >( dMicroGradientYielddDrivingStress.data( ) );
            auto map_dMicroGradientYielddCohesion      = getFixedSizeMatrixMap< floatType,       3,       3 >( dMicroGradientYielddCohesion.data( ) );
            auto map_dMicroGradientYielddPrecedingF    = getFixedSizeMatrixMap< floatType,       3, sot_dim >( dMicroGradientYielddPrecedingF.data( ) );

            auto map_dPrecedingFdF                         = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPrecedingFdF->data( ) );
            auto map_dMacroDrivingStressdStress            = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMacroDrivingStressdStress->data( ) );
            auto map_dMicroDrivingStressdStress            = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMicroDrivingStressdStress->data( ) );
            auto map_dMicroGradientDrivingStressdStress    = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dMicroGradientDrivingStressdStress->data( ) );
            auto map_dMacroDrivingStressdF                 = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMacroDrivingStressdF->data( ) );
            auto map_dMicroDrivingStressdF                 = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMicroDrivingStressdF->data( ) );
            auto map_dMicroGradientDrivingStressdF         = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dMicroGradientDrivingStressdF->data( ) );
            auto map_dMicroGradientDrivingStressdChi       = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dMicroGradientDrivingStressdChi->data( ) );

            auto map_dPrecedingFdFn                        = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dPrecedingFdFn->data( ),                        ( num_configs - 1 ) * sot_dim );
            auto map_dMicroGradientCohesiondStateVariables = getDynamicColumnSizeMatrixMap< floatType,       3 >( dMicroGradientCohesiondStateVariables->data( ), num_isvs                      );
            auto map_dMacroDrivingStressdFn                = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dMacroDrivingStressdFn->data( ),                ( num_configs - 1 ) * sot_dim    );
            auto map_dMicroDrivingStressdFn                = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dMicroDrivingStressdFn->data( ),                ( num_configs - 1 ) * sot_dim    );
            auto map_dMicroGradientDrivingStressdFn        = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dMicroGradientDrivingStressdFn->data( ),        ( num_configs - 1 ) * sot_dim    );
            auto map_dMicroGradientDrivingStressdChin      = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dMicroGradientDrivingStressdChin->data( ),      ( num_configs - 1 ) * sot_dim    );

            auto map_dMacroYielddStress                 = dMacroYielddStress.zeroMap< floatType, 1, sot_dim >( );
            auto map_dMacroYielddF                      = dMacroYielddF.zeroMap< floatType, 1, sot_dim >( );
            auto map_dMacroYielddFn                     = dMacroYielddFn.zeroMap< floatType, 1 >( ( num_configs - 1 ) * sot_dim );
            auto map_dMicroYielddStress                 = dMicroYielddStress.zeroMap< floatType, 1, sot_dim >( );
            auto map_dMicroYielddF                      = dMicroYielddF.zeroMap< floatType, 1, sot_dim >( );
            auto map_dMicroYielddFn                     = dMicroYielddFn.zeroMap< floatType, 1 >( ( num_configs - 1 ) * sot_dim );
            auto map_dMicroGradientYielddStress         = dMicroGradientYielddStress.zeroMap< floatType, 3, tot_dim >( );
            auto map_dMicroGradientYielddStateVariables = dMicroGradientYielddStateVariables.zeroMap< floatType, 3 >( ( num_configs - 1 ) * num_isvs );
            auto map_dMicroGradientYielddF              = dMicroGradientYielddF.zeroMap< floatType, 3, sot_dim >( );
            auto map_dMicroGradientYielddFn             = dMicroGradientYielddFn.zeroMap< floatType, 3 >( ( num_configs - 1 ) * sot_dim );
            auto map_dMicroGradientYielddChi            = dMicroGradientYielddChi.zeroMap< floatType, 3, sot_dim >( );
            auto map_dMicroGradientYielddChin           = dMicroGradientYielddChin.zeroMap< floatType, 3 >( ( num_configs - 1 ) * sot_dim );

            map_dMacroYielddStress = ( map_dMacroYielddDrivingStress * map_dMacroDrivingStressdStress ).eval( );

            *dMacroYielddStateVariables.value = dMacroYielddCohesion * ( *dMacroCohesiondStateVariables );

            map_dMacroYielddF  = ( map_dMacroYielddDrivingStress * map_dMacroDrivingStressdF ).eval( );
            map_dMacroYielddF += ( map_dMacroYielddPrecedingF * map_dPrecedingFdF ).eval( );

            map_dMacroYielddFn  = ( map_dMacroYielddDrivingStress * map_dMacroDrivingStressdFn ).eval( );
            map_dMacroYielddFn += ( map_dMacroYielddPrecedingF * map_dPrecedingFdFn ).eval( );

            map_dMicroYielddStress = ( map_dMicroYielddDrivingStress * map_dMicroDrivingStressdStress ).eval( );

            *dMicroYielddStateVariables.value = dMicroYielddCohesion * ( *dMicroCohesiondStateVariables );

            map_dMicroYielddF  = ( map_dMicroYielddDrivingStress * map_dMicroDrivingStressdF ).eval( );
            map_dMicroYielddF += ( map_dMicroYielddPrecedingF * map_dPrecedingFdF ).eval( );

            map_dMicroYielddFn  = ( map_dMicroYielddDrivingStress * map_dMicroDrivingStressdFn ).eval( );
            map_dMicroYielddFn += ( map_dMicroYielddPrecedingF * map_dPrecedingFdFn ).eval( );

            map_dMicroGradientYielddStress = ( map_dMicroGradientYielddDrivingStress * map_dMicroGradientDrivingStressdStress ).eval( );

            map_dMicroGradientYielddStateVariables = ( map_dMicroGradientYielddCohesion * map_dMicroGradientCohesiondStateVariables ).eval( );

            map_dMicroGradientYielddF  = ( map_dMicroGradientYielddDrivingStress * map_dMicroGradientDrivingStressdF ).eval( );
            map_dMicroGradientYielddF += ( map_dMicroGradientYielddPrecedingF * map_dPrecedingFdF ).eval( );

            map_dMicroGradientYielddFn  = ( map_dMicroGradientYielddDrivingStress * map_dMicroGradientDrivingStressdFn ).eval( );
            map_dMicroGradientYielddFn += ( map_dMicroGradientYielddPrecedingF * map_dPrecedingFdFn ).eval( );

            map_dMicroGradientYielddChi = ( map_dMicroGradientYielddDrivingStress * map_dMicroGradientDrivingStressdChi ).eval( );

            map_dMicroGradientYielddChin = ( map_dMicroGradientYielddDrivingStress * map_dMicroGradientDrivingStressdChin ).eval( );
//            thirdOrderTensor dMicroGradientYielddF = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdF, 3, tot_dim, tot_dim, sot_dim )
//                                                   + tardigradeVectorTools::matrixMultiply( dMicroGradientYielddPrecedingF, *dPrecedingFdF, 3, sot_dim, sot_dim, sot_dim );
//
//            floatVector dMicroGradientYielddFn = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdFn, 3, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
//                                               + tardigradeVectorTools::matrixMultiply( dMicroGradientYielddPrecedingF, *dPrecedingFdFn, 3, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );
//
//            thirdOrderTensor dMicroGradientYielddChi = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChi, 3, tot_dim, tot_dim, sot_dim );
//
//            floatVector dMicroGradientYielddChin = tardigradeVectorTools::matrixMultiply( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChin, 3, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            if ( isPrevious ){

//                set_previousdMicroGradientYielddF( dMicroGradientYielddF );
//
//                set_previousdMicroGradientYielddFn( dMicroGradientYielddFn );
//
//                set_previousdMicroGradientYielddChi( dMicroGradientYielddChi );
//
//                set_previousdMicroGradientYielddChin( dMicroGradientYielddChin );

            }
            else{

//                set_dMicroGradientYielddF( dMicroGradientYielddF );
//
//                set_dMicroGradientYielddFn( dMicroGradientYielddFn );
//
//                set_dMicroGradientYielddChi( dMicroGradientYielddChi );
//
//                set_dMicroGradientYielddChin( dMicroGradientYielddChin );

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

            const fourthOrderTensor *dF1dF;

            const floatVector *dF1dFn;

            setDataStorageBase< fourthOrderTensor > dPrecedingFdF;

            setDataStorageBase< floatVector > dPrecedingFdFn;

            if ( isPrevious ){

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF  = hydra->get_previousdF1dF( );

                dF1dFn = hydra->get_previousdF1dFn( );

                dPrecedingFdF  = get_setDataStorage_previousdPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_setDataStorage_previousdPrecedingDeformationGradientdFn( );

            }
            else{

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF  = hydra->get_dF1dF( );

                dF1dFn = hydra->get_dF1dFn( );

                dPrecedingFdF  = get_setDataStorage_dPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_setDataStorage_dPrecedingDeformationGradientdFn( );

            }

            // Construct the derivatives of the preceding F

            dPrecedingFdF.zero( sot_dim * sot_dim );

            dPrecedingFdFn.zero( sot_dim * ( num_configs - 1 ) * sot_dim );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dPrecedingFdF.value )[ sot_dim * i + j ] += dPrecedingFdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dF )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: See if this can be sped up by accessing the cache better

                    ( *dPrecedingFdFn.value )[ ( num_configs - 1 ) * sot_dim * i + j ] = dPrecedingFdSubFs[ num_configs * sot_dim * i + sot_dim + j ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dPrecedingFdFn.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dPrecedingFdSubFs[ num_configs * sot_dim * i + k ] * ( *dF1dFn )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

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

            const fourthOrderTensor *dChi1dChi;

            const floatVector *dChi1dChin;

            setDataStorageBase< fourthOrderTensor > dPrecedingChidChi;

            setDataStorageBase< floatVector >       dPrecedingChidChin;

            if ( isPrevious ){

                set_previousPrecedingMicroDeformation( hydra->getPreviousPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingChidSubChis = hydra->getPreviousPrecedingMicroConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dChi1dChi  = hydra->get_previousdChi1dChi( );

                dChi1dChin = hydra->get_previousdChi1dChin( );

                dPrecedingChidChi  = get_setDataStorage_previousdPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_setDataStorage_previousdPrecedingMicroDeformationdChin( );

            }
            else{

                set_precedingMicroDeformation( hydra->getPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingChidSubChis = hydra->getPrecedingMicroConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dChi1dChi  = hydra->get_dChi1dChi( );

                dChi1dChin = hydra->get_dChi1dChin( );

                dPrecedingChidChi  = get_setDataStorage_dPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_setDataStorage_dPrecedingMicroDeformationdChin( );

            }

            // Construct the derivatives of the preceding F

            dPrecedingChidChi.zero( sot_dim * sot_dim );
            dPrecedingChidChin.zero( sot_dim * sot_dim * ( num_configs - 1 ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dPrecedingChidChi.value )[ sot_dim * i + j ] += dPrecedingChidSubChis[ num_configs * sot_dim * i + k ] * ( *dChi1dChi )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: See if this can be sped up by accessing the cache better

                    ( *dPrecedingChidChin.value )[ ( num_configs - 1 ) * sot_dim * i + j ] = dPrecedingChidSubChis[ num_configs * sot_dim * i + sot_dim + j ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dPrecedingChidChin.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dPrecedingChidSubChis[ num_configs * sot_dim * i + k ] * ( *dChi1dChin )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

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

                set_previousPrecedingGradientMicroDeformation( thirdOrderTensor( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                                                 hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression

            }
            else{

                set_precedingGradientMicroDeformation( thirdOrderTensor( hydra->get_gradientMicroConfigurations( )->begin( ),
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

                set_previousPrecedingGradientMicroDeformation( thirdOrderTensor( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                                                 hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim ) ); //TODO: Generalize this expression

                set_previousdPrecedingGradientMicroDeformationdFn(       *hydra->get_previousdGradChi1dFn( ) );

                set_previousdPrecedingGradientMicroDeformationdChi(      *hydra->get_previousdGradChi1dChi( ) );

                set_previousdPrecedingGradientMicroDeformationdChin(     *hydra->get_previousdGradChi1dChin( ) );

                set_previousdPrecedingGradientMicroDeformationdGradChi(  *hydra->get_previousdGradChi1dGradChi( ) );

                set_previousdPrecedingGradientMicroDeformationdGradChin( *hydra->get_previousdGradChi1dGradChin( ) );

            }
            else{

                set_precedingGradientMicroDeformation( thirdOrderTensor( hydra->get_gradientMicroConfigurations( )->begin( ),
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

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const secondOrderTensor *precedingDeformationGradient;

            const secondOrderTensor *precedingMicroDeformation;

            const thirdOrderTensor  *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const secondOrderTensor *dMacroFlowdDrivingStress;

            const secondOrderTensor *dMicroFlowdDrivingStress;

            const fourthOrderTensor *dMicroGradientFlowdDrivingStress;

            setDataStorageBase< secondOrderTensor > macroVelocityGradient;

            setDataStorageBase< secondOrderTensor > microVelocityGradient;

            setDataStorageBase< thirdOrderTensor >  gradientMicroVelocityGradient;

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                precedingMicroDeformation = get_previousPrecedingMicroDeformation( );

                precedingGradientMicroDeformation = get_previousPrecedingGradientMicroDeformation( );

                plasticMultipliers = get_previousPlasticMultipliers( );

                dMacroFlowdDrivingStress = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_previousdMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_previousdMicroGradientFlowdDrivingStress( );

                macroVelocityGradient         = get_setDataStorage_previousPlasticMacroVelocityGradient( );

                microVelocityGradient         = get_setDataStorage_previousPlasticMicroVelocityGradient( );

                gradientMicroVelocityGradient = get_setDataStorage_previousPlasticGradientMicroVelocityGradient( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

                precedingMicroDeformation = get_precedingMicroDeformation( );

                precedingGradientMicroDeformation = get_precedingGradientMicroDeformation( );

                plasticMultipliers = get_plasticMultipliers( );

                dMacroFlowdDrivingStress = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_dMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_dMicroGradientFlowdDrivingStress( );

                macroVelocityGradient         = get_setDataStorage_plasticMacroVelocityGradient( );

                microVelocityGradient         = get_setDataStorage_plasticMicroVelocityGradient( );

                gradientMicroVelocityGradient = get_setDataStorage_plasticGradientMicroVelocityGradient( );

            }

            // Form the preceding RCG and its inverse
            secondOrderTensor precedingRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG ) );

            secondOrderTensor inversePrecedingRCG = precedingRCG;
            auto map_inversePrecedingRCG = getFixedSizeMatrixMap< floatType, dim, dim >( inversePrecedingRCG.data( ) );
            map_inversePrecedingRCG = map_inversePrecedingRCG.inverse( ).eval( );

            // Form the precedingPsi and its inverse
            secondOrderTensor precedingPsi( sot_dim, 0 );
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( *precedingDeformationGradient, *precedingMicroDeformation, precedingPsi ) );

            secondOrderTensor inversePrecedingPsi = precedingPsi;
            auto map_inversePrecedingPsi = getFixedSizeMatrixMap< floatType, dim, dim >( inversePrecedingPsi.data( ) );
            map_inversePrecedingPsi = map_inversePrecedingPsi.inverse( ).eval( );

            // Form the preceding micro RCG and its inverse
            secondOrderTensor precedingMicroRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG ) );

            // Form Gamma
            thirdOrderTensor precedingGamma;
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( *precedingDeformationGradient, *precedingGradientMicroDeformation, precedingGamma ) );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( ( *plasticMultipliers )[ 0 ], ( *plasticMultipliers )[ 1 ],
                                                     inversePrecedingRCG, *dMacroFlowdDrivingStress, *dMicroFlowdDrivingStress,
                                                     *macroVelocityGradient.value )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( ( *plasticMultipliers )[ 1 ], precedingMicroRCG, precedingPsi, inversePrecedingPsi,
                                                     *dMicroFlowdDrivingStress, *microVelocityGradient.value );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroGradientVelocityGradient( floatVector( plasticMultipliers->begin( ) + 2, plasticMultipliers->end( ) ),
                                                             precedingPsi, inversePrecedingPsi, precedingGamma,
                                                             *dMicroGradientFlowdDrivingStress,
                                                             *microVelocityGradient.value, *gradientMicroVelocityGradient.value );
            )

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

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            constexpr unsigned int fiot_dim = fot_dim * dim;

            constexpr unsigned int siot_dim = fiot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const secondOrderTensor *precedingDeformationGradient;

            const fourthOrderTensor *dPrecedingFdF;

            const floatVector *dPrecedingFdFn;

            const fourthOrderTensor *dPrecedingChidChi;

            const floatVector *dPrecedingChidChin;

            const floatVector *dPrecedingGradChidFn;

            const fifthOrderTensor *dPrecedingGradChidChi;

            const floatVector *dPrecedingGradChidChin;

            const sixthOrderTensor *dPrecedingGradChidGradChi;

            const floatVector *dPrecedingGradChidGradChin;

            const secondOrderTensor *precedingMicroDeformation;

            const thirdOrderTensor *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const secondOrderTensor *dMacroFlowdDrivingStress;

            const secondOrderTensor *dMicroFlowdDrivingStress;

            const fourthOrderTensor *dMicroGradientFlowdDrivingStress;

            const seventhOrderTensor *d2MacroFlowdDrivingStressdStress;

            const fourthOrderTensor *d2MacroFlowdDrivingStressdF;

            const floatVector *d2MacroFlowdDrivingStressdFn;

            const fourthOrderTensor *d2MicroFlowdDrivingStressdF;

            const floatVector *d2MicroFlowdDrivingStressdFn;

            const fourthOrderTensor *d2MicroFlowdDrivingStressdStress;

            const seventhOrderTensor *d2MicroGradientFlowdDrivingStressdStress;

            const sixthOrderTensor *d2MicroGradientFlowdDrivingStressdF;

            const floatVector *d2MicroGradientFlowdDrivingStressdFn;

            const sixthOrderTensor *d2MicroGradientFlowdDrivingStressdChi;

            const floatVector *d2MicroGradientFlowdDrivingStressdChin;

            setDataStorageBase< secondOrderTensor > macroVelocityGradient;

            setDataStorageBase< secondOrderTensor > microVelocityGradient;

            setDataStorageBase< thirdOrderTensor >  gradientMicroVelocityGradient;

            setDataStorageBase< secondOrderTensor > dPlasticMacroLdMacroStress;

            setDataStorageBase< secondOrderTensor > dPlasticMacroLdMicroStress;

            setDataStorageBase< secondOrderTensor > dPlasticMicroLdMicroStress;

            setDataStorageBase< fifthOrderTensor >  dPlasticGradientMicroLdMicroStress;

            setDataStorageBase< sixthOrderTensor >  dPlasticGradientMicroLdHigherOrderStress;

            setDataStorageBase< fourthOrderTensor > dPlasticMacroLdF;

            setDataStorageBase< floatVector > dPlasticMacroLdFn;

            setDataStorageBase< fourthOrderTensor > dPlasticMicroLdF;

            setDataStorageBase< floatVector > dPlasticMicroLdFn;

            setDataStorageBase< fifthOrderTensor > dPlasticGradientMicroLdF;

            setDataStorageBase< floatVector > dPlasticGradientMicroLdFn;

            setDataStorageBase< fourthOrderTensor > dPlasticMicroLdChi;

            setDataStorageBase< floatVector > dPlasticMicroLdChin;

            setDataStorageBase< fifthOrderTensor > dPlasticGradientMicroLdChi;

            setDataStorageBase< floatVector > dPlasticGradientMicroLdChin;

            setDataStorageBase< sixthOrderTensor > dPlasticGradientMicroLdGradChi;

            setDataStorageBase< floatVector > dPlasticGradientMicroLdGradChin;

            setDataStorageBase< floatVector > dPlasticMacroLdISVs;

            setDataStorageBase< floatVector > dPlasticMicroLdISVs;

            setDataStorageBase< floatVector > dPlasticGradientMicroLdISVs;

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

                d2MicroGradientFlowdDrivingStressdStress = get_previousd2MicroGradientFlowdDrivingStressdStress( );

                d2MicroGradientFlowdDrivingStressdF      = get_previousd2MicroGradientFlowdDrivingStressdF( );

                d2MicroGradientFlowdDrivingStressdFn     = get_previousd2MicroGradientFlowdDrivingStressdFn( );

                d2MicroGradientFlowdDrivingStressdChi    = get_previousd2MicroGradientFlowdDrivingStressdChi( );

                d2MicroGradientFlowdDrivingStressdChin   = get_previousd2MicroGradientFlowdDrivingStressdChin( );

                dMacroFlowdDrivingStress         = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_previousdMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_previousdMicroGradientFlowdDrivingStress( );

                macroVelocityGradient         = get_setDataStorage_previousPlasticMacroVelocityGradient( );

                microVelocityGradient         = get_setDataStorage_previousPlasticMicroVelocityGradient( );

                gradientMicroVelocityGradient = get_setDataStorage_previousPlasticGradientMicroVelocityGradient( );

                dPlasticMacroLdMacroStress    = get_setDataStorage_previousdPlasticMacroVelocityGradientdMacroStress( );

                dPlasticMacroLdMicroStress    = get_setDataStorage_previousdPlasticMacroVelocityGradientdMicroStress( );

                dPlasticMicroLdMicroStress    = get_setDataStorage_previousdPlasticMicroVelocityGradientdMicroStress( );

                dPlasticGradientMicroLdMicroStress = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdMicroStress( );

                dPlasticGradientMicroLdHigherOrderStress = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( );

                dPlasticMacroLdF  = get_setDataStorage_previousdPlasticMacroVelocityGradientdF( );

                dPlasticMacroLdFn = get_setDataStorage_previousdPlasticMacroVelocityGradientdFn( );

                dPlasticMicroLdF  = get_setDataStorage_previousdPlasticMicroVelocityGradientdF( );

                dPlasticMicroLdFn = get_setDataStorage_previousdPlasticMicroVelocityGradientdFn( );

                dPlasticGradientMicroLdF  = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdF( );

                dPlasticGradientMicroLdFn = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdFn( );

                dPlasticMicroLdChi  = get_setDataStorage_previousdPlasticMicroVelocityGradientdChi( );

                dPlasticMicroLdChin = get_setDataStorage_previousdPlasticMicroVelocityGradientdChin( );

                dPlasticGradientMicroLdChi  = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdChi( );

                dPlasticGradientMicroLdChin = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdChin( );

                dPlasticGradientMicroLdGradChi  = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdGradChi( );

                dPlasticGradientMicroLdGradChin = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdGradChin( );

                dPlasticMacroLdISVs = get_setDataStorage_previousdPlasticMacroVelocityGradientdStateVariables( );

                dPlasticMicroLdISVs = get_setDataStorage_previousdPlasticMicroVelocityGradientdStateVariables( );

                dPlasticGradientMicroLdISVs = get_setDataStorage_previousdPlasticGradientMicroVelocityGradientdStateVariables( );

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

                d2MicroGradientFlowdDrivingStressdStress = get_d2MicroGradientFlowdDrivingStressdStress( );

                d2MicroGradientFlowdDrivingStressdF      = get_d2MicroGradientFlowdDrivingStressdF( );

                d2MicroGradientFlowdDrivingStressdFn     = get_d2MicroGradientFlowdDrivingStressdFn( );

                d2MicroGradientFlowdDrivingStressdChi    = get_d2MicroGradientFlowdDrivingStressdChi( );

                d2MicroGradientFlowdDrivingStressdChin   = get_d2MicroGradientFlowdDrivingStressdChin( );

                dMacroFlowdDrivingStress         = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_dMicroFlowdDrivingStress( );

                dMicroGradientFlowdDrivingStress = get_dMicroGradientFlowdDrivingStress( );

                macroVelocityGradient         = get_setDataStorage_plasticMacroVelocityGradient( );

                microVelocityGradient         = get_setDataStorage_plasticMicroVelocityGradient( );

                gradientMicroVelocityGradient = get_setDataStorage_plasticGradientMicroVelocityGradient( );

                dPlasticMacroLdMacroStress    = get_setDataStorage_dPlasticMacroVelocityGradientdMacroStress( );

                dPlasticMacroLdMicroStress    = get_setDataStorage_dPlasticMacroVelocityGradientdMicroStress( );

                dPlasticMicroLdMicroStress    = get_setDataStorage_dPlasticMicroVelocityGradientdMicroStress( );

                dPlasticGradientMicroLdMicroStress = get_setDataStorage_dPlasticGradientMicroVelocityGradientdMicroStress( );

                dPlasticGradientMicroLdHigherOrderStress = get_setDataStorage_dPlasticGradientMicroVelocityGradientdHigherOrderStress( );

                dPlasticMacroLdF  = get_setDataStorage_dPlasticMacroVelocityGradientdF( );

                dPlasticMacroLdFn = get_setDataStorage_dPlasticMacroVelocityGradientdFn( );

                dPlasticMicroLdF  = get_setDataStorage_dPlasticMicroVelocityGradientdF( );

                dPlasticMicroLdFn = get_setDataStorage_dPlasticMicroVelocityGradientdFn( );

                dPlasticGradientMicroLdF  = get_setDataStorage_dPlasticGradientMicroVelocityGradientdF( );

                dPlasticGradientMicroLdFn = get_setDataStorage_dPlasticGradientMicroVelocityGradientdFn( );

                dPlasticMicroLdChi  = get_setDataStorage_dPlasticMicroVelocityGradientdChi( );

                dPlasticMicroLdChin = get_setDataStorage_dPlasticMicroVelocityGradientdChin( );

                dPlasticGradientMicroLdChi  = get_setDataStorage_dPlasticGradientMicroVelocityGradientdChi( );

                dPlasticGradientMicroLdChin = get_setDataStorage_dPlasticGradientMicroVelocityGradientdChin( );

                dPlasticGradientMicroLdGradChi  = get_setDataStorage_dPlasticGradientMicroVelocityGradientdGradChi( );

                dPlasticGradientMicroLdGradChin = get_setDataStorage_dPlasticGradientMicroVelocityGradientdGradChin( );

                dPlasticMacroLdISVs = get_setDataStorage_dPlasticMacroVelocityGradientdStateVariables( );

                dPlasticMicroLdISVs = get_setDataStorage_dPlasticMicroVelocityGradientdStateVariables( );

                dPlasticGradientMicroLdISVs = get_setDataStorage_dPlasticGradientMicroVelocityGradientdStateVariables( );

            }

            // Form the preceding RCG and its inverse
            secondOrderTensor precedingRCG;

            secondOrderTensor dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG, dRCGdPrecedingF ) );

            secondOrderTensor inversePrecedingRCG = precedingRCG;
            auto map_inversePrecedingRCG = getFixedSizeMatrixMap< floatType, dim, dim >( inversePrecedingRCG.data( ) );
            map_inversePrecedingRCG = map_inversePrecedingRCG.inverse( ).eval( );

            auto map_dRCGdPrecedingF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dRCGdPrecedingF.data( ) );
            auto map_dPrecedingFdF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPrecedingFdF->data( ) );
            auto map_dPrecedingFdFn = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dPrecedingFdFn->data( ), ( num_configs - 1 ) * sot_dim );

            fourthOrderTensor dRCGdF( fot_dim, 0 );
            auto map_dRCGdF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dRCGdF.data( ) );

            floatVector dRCGdFn( fot_dim * ( num_configs - 1 ), 0 );
            auto map_dRCGdFn = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dRCGdFn.data( ), ( num_configs - 1 ) * sot_dim );

            map_dRCGdF = ( map_dRCGdPrecedingF * map_dPrecedingFdF ).eval( );

            map_dRCGdFn = ( map_dRCGdPrecedingF * map_dPrecedingFdFn ).eval( );

            // Form the precedingPsi and its inverse
            secondOrderTensor precedingPsi;

            fourthOrderTensor dPsidPrecedingF, dPsidPrecedingChi;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( *precedingDeformationGradient, *precedingMicroDeformation, precedingPsi, dPsidPrecedingF, dPsidPrecedingChi ) );

            auto map_dPsidPrecedingF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPsidPrecedingF.data( ) );
            auto map_dPsidPrecedingChi = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPsidPrecedingChi.data( ) );
            auto map_dPrecedingChidChi = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPrecedingChidChi->data( ) );
            auto map_dPrecedingChidChin = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dPrecedingChidChin->data( ), ( num_configs - 1 ) * sot_dim );

            fourthOrderTensor dPsidF( fot_dim, 0 );
            auto map_dPsidF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPsidF.data( ) );

            fourthOrderTensor dPsidFn( fot_dim * ( num_configs - 1 ), 0 );
            auto map_dPsidFn = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dPsidFn.data( ), ( num_configs - 1 ) * sot_dim );

            fourthOrderTensor dPsidChi( fot_dim, 0 );
            auto map_dPsidChi = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPsidChi.data( ) );

            fourthOrderTensor dPsidChin( fot_dim * ( num_configs - 1 ), 0 );
            auto map_dPsidChin = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dPsidChin.data( ), ( num_configs - 1 ) * sot_dim );

            map_dPsidF = ( map_dPsidPrecedingF * map_dPrecedingFdF ).eval( );

            map_dPsidFn = ( map_dPsidPrecedingF * map_dPrecedingFdFn ).eval( );

            map_dPsidChi = ( map_dPsidPrecedingChi * map_dPrecedingChidChi ).eval( );

            map_dPsidChin = ( map_dPsidPrecedingChi * map_dPrecedingChidChin ).eval( );

            secondOrderTensor inversePrecedingPsi = precedingPsi;
            auto map_inversePrecedingPsi = getFixedSizeMatrixMap< floatType, dim, dim >( inversePrecedingPsi.data( ) );
            map_inversePrecedingPsi = map_inversePrecedingPsi.inverse( ).eval( );

            // Form the preceding micro RCG and its inverse
            secondOrderTensor precedingMicroRCG;

            fourthOrderTensor dMicroRCGdPrecedingChi;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG, dMicroRCGdPrecedingChi ) );

            auto map_dMicroRCGdPrecedingChi = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMicroRCGdPrecedingChi.data( ) );

            fourthOrderTensor dMicroRCGdChi( fot_dim, 0 );
            auto map_dMicroRCGdChi = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dMicroRCGdChi.data( ) );

            fourthOrderTensor dMicroRCGdChin( fot_dim * ( num_configs - 1 ), 0 );
            auto map_dMicroRCGdChin = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dMicroRCGdChin.data( ), ( num_configs - 1 ) * sot_dim );

            map_dMicroRCGdChi = ( map_dMicroRCGdPrecedingChi * map_dPrecedingChidChi ).eval( );

            map_dMicroRCGdChin = ( map_dMicroRCGdPrecedingChi * map_dPrecedingChidChin ).eval( );

            // Form Gamma
            thirdOrderTensor precedingGamma;

            fifthOrderTensor dPrecedingGammadPrecedingF;

            sixthOrderTensor dPrecedingGammadPrecedingGradChi;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( *precedingDeformationGradient, *precedingGradientMicroDeformation, precedingGamma, dPrecedingGammadPrecedingF, dPrecedingGammadPrecedingGradChi ) );

            auto map_dPrecedingGammadPrecedingF = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPrecedingGammadPrecedingF.data( ) );
            auto map_dPrecedingGammadPrecedingGradChi = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dPrecedingGammadPrecedingGradChi.data( ) );
            auto map_dPrecedingGradChidFn = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGradChidFn->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_dPrecedingGradChidChi  = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPrecedingGradChidChi->data( ) );
            auto map_dPrecedingGradChidChin = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGradChidChin->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_dPrecedingGradChidGradChi  = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dPrecedingGradChidGradChi->data( ) );
            auto map_dPrecedingGradChidGradChin = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGradChidGradChin->data( ), ( num_configs - 1 ) * tot_dim );

            fifthOrderTensor dPrecedingGammadF( fiot_dim, 0 );
            auto map_dPrecedingGammadF = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPrecedingGammadF.data( ) );

            floatVector dPrecedingGammadFn( fiot_dim * ( num_configs - 1 ), 0 );
            auto map_dPrecedingGammadFn = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGammadFn.data( ), ( num_configs - 1 ) * sot_dim );

            fifthOrderTensor dPrecedingGammadChi( fiot_dim, 0 );
            auto map_dPrecedingGammadChi = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPrecedingGammadChi.data( ) );

            floatVector dPrecedingGammadChin( fiot_dim * ( num_configs - 1 ), 0 );
            auto map_dPrecedingGammadChin = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGammadChin.data( ), ( num_configs - 1 ) * sot_dim );

            sixthOrderTensor dPrecedingGammadGradChi( siot_dim, 0 );
            auto map_dPrecedingGammadGradChi = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dPrecedingGammadGradChi.data( ) );

            floatVector dPrecedingGammadGradChin( siot_dim * ( num_configs - 1 ), 0 );
            auto map_dPrecedingGammadGradChin = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( dPrecedingGammadGradChin.data( ), ( num_configs - 1 ) * tot_dim );

            map_dPrecedingGammadF = ( map_dPrecedingGammadPrecedingF * map_dPrecedingFdF ).eval( );

            map_dPrecedingGammadFn  = ( map_dPrecedingGammadPrecedingF * map_dPrecedingFdFn ).eval( );
            map_dPrecedingGammadFn += ( map_dPrecedingGammadPrecedingGradChi * map_dPrecedingGradChidFn ).eval( );

            map_dPrecedingGammadChi = ( map_dPrecedingGammadPrecedingGradChi * map_dPrecedingGradChidChi ).eval( );

            map_dPrecedingGammadChin = ( map_dPrecedingGammadPrecedingGradChi * map_dPrecedingGradChidChin ).eval( );

            map_dPrecedingGammadGradChi = ( map_dPrecedingGammadPrecedingGradChi * map_dPrecedingGradChidGradChi ).eval( );

            map_dPrecedingGammadGradChin = ( map_dPrecedingGammadPrecedingGradChi * map_dPrecedingGradChidGradChin ).eval( );

            secondOrderTensor dPlasticMacroLdMacroGamma;

            secondOrderTensor dPlasticMacroLdMicroGamma;

            fourthOrderTensor dPlasticMacroLdPrecedingRCG;

            fourthOrderTensor dPlasticMacroLdMacroFlowDirection;

            fourthOrderTensor dPlasticMacroLdMicroFlowDirection;

            secondOrderTensor dPlasticMicroLdMicroGamma;

            fourthOrderTensor dPlasticMicroLdPrecedingMicroRCG;

            fourthOrderTensor dPlasticMicroLdPrecedingPsi;

            fourthOrderTensor dPlasticMicroLdMicroFlowDirection;

            fourthOrderTensor dPlasticGradientMicroLdMicroGradientGamma;

            fifthOrderTensor  dPlasticGradientMicroLdPlasticMicroL;

            fifthOrderTensor  dPlasticGradientMicroLdPrecedingPsi;

            sixthOrderTensor  dPlasticGradientMicroLdPrecedingGamma;

            sixthOrderTensor  dPlasticGradientMicroLdMicroGradientFlowDirection;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( ( *plasticMultipliers )[ 0 ], ( *plasticMultipliers )[ 1 ],
                                                     inversePrecedingRCG, *dMacroFlowdDrivingStress, *dMicroFlowdDrivingStress,
                                                     *macroVelocityGradient.value, dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma,
                                                     dPlasticMacroLdPrecedingRCG, dPlasticMacroLdMacroFlowDirection, dPlasticMacroLdMicroFlowDirection )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( ( *plasticMultipliers )[ 1 ], precedingMicroRCG, precedingPsi, inversePrecedingPsi,
                                                     *dMicroFlowdDrivingStress, *microVelocityGradient.value, dPlasticMicroLdMicroGamma,
                                                     dPlasticMicroLdPrecedingMicroRCG, dPlasticMicroLdPrecedingPsi, dPlasticMicroLdMicroFlowDirection );
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroGradientVelocityGradient( floatVector( plasticMultipliers->begin( ) + 2, plasticMultipliers->end( ) ),
                                                             precedingPsi, inversePrecedingPsi, precedingGamma,
                                                             *dMicroGradientFlowdDrivingStress,
                                                             *microVelocityGradient.value, *gradientMicroVelocityGradient.value,
                                                             dPlasticGradientMicroLdMicroGradientGamma,
                                                             dPlasticGradientMicroLdPlasticMicroL,
                                                             dPlasticGradientMicroLdPrecedingPsi,
                                                             dPlasticGradientMicroLdPrecedingGamma,
                                                             dPlasticGradientMicroLdMicroGradientFlowDirection );
            )

            auto map_dPlasticMacroLdPrecedingRCG                       = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMacroLdPrecedingRCG.data( ) );
            auto map_dPlasticMacroLdMacroFlowDirection                 = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMacroLdMacroFlowDirection.data( ) );
            auto map_dPlasticMacroLdMicroFlowDirection                 = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMacroLdMicroFlowDirection.data( ) );

            auto map_dPlasticMicroLdPrecedingMicroRCG                  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMicroLdPrecedingMicroRCG.data( ) );
            auto map_dPlasticMicroLdPrecedingPsi                       = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMicroLdPrecedingPsi.data( ) );
            auto map_dPlasticMicroLdMicroFlowDirection                 = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMicroLdMicroFlowDirection.data( ) );

            auto map_dPlasticGradientMicroLdPlasticMicroL              = getFixedSizeMatrixMap< floatType, tot_dim,     sot_dim >( dPlasticGradientMicroLdPlasticMicroL.data( ) );
            auto map_dPlasticGradientMicroLdPrecedingPsi               = getFixedSizeMatrixMap< floatType, tot_dim,     sot_dim >( dPlasticGradientMicroLdPrecedingPsi.data( ) );
            auto map_dPlasticGradientMicroLdPrecedingGamma             = getFixedSizeMatrixMap< floatType, tot_dim,     tot_dim >( dPlasticGradientMicroLdPrecedingGamma.data( ) );
            auto map_dPlasticGradientMicroLdMicroGradientFlowDirection = getFixedSizeMatrixMap< floatType, tot_dim, 3 * tot_dim >( dPlasticGradientMicroLdMicroGradientFlowDirection.data( ) );

            auto map_d2MacroFlowdDrivingStressdStress         = getFixedSizeMatrixMap< floatType,     sot_dim, sot_dim >( d2MacroFlowdDrivingStressdStress->data( ) );
            auto map_d2MicroFlowdDrivingStressdStress         = getFixedSizeMatrixMap< floatType,     sot_dim, sot_dim >( d2MicroFlowdDrivingStressdStress->data( ) );
            auto map_d2MicroGradientFlowdDrivingStressdStress = getFixedSizeMatrixMap< floatType, 3 * tot_dim, tot_dim >( d2MicroGradientFlowdDrivingStressdStress->data( ) );
            auto map_d2MacroFlowdDrivingStressdF              = getFixedSizeMatrixMap< floatType,     sot_dim, sot_dim >( d2MacroFlowdDrivingStressdF->data( ) );
            auto map_d2MicroFlowdDrivingStressdF              = getFixedSizeMatrixMap< floatType,     sot_dim, sot_dim >( d2MicroFlowdDrivingStressdF->data( ) );
            auto map_d2MicroGradientFlowdDrivingStressdF      = getFixedSizeMatrixMap< floatType, 3 * tot_dim, sot_dim >( d2MicroGradientFlowdDrivingStressdF->data( ) );
            auto map_d2MicroGradientFlowdDrivingStressdChi    = getFixedSizeMatrixMap< floatType, 3 * tot_dim, sot_dim >( d2MicroGradientFlowdDrivingStressdChi->data( ) );
            auto map_d2MacroFlowdDrivingStressdFn             = getDynamicColumnSizeMatrixMap< floatType,     sot_dim >( d2MacroFlowdDrivingStressdFn->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_d2MicroFlowdDrivingStressdFn             = getDynamicColumnSizeMatrixMap< floatType,     sot_dim >( d2MicroFlowdDrivingStressdFn->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_d2MicroGradientFlowdDrivingStressdFn     = getDynamicColumnSizeMatrixMap< floatType, 3 * tot_dim >( d2MicroGradientFlowdDrivingStressdFn->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_d2MicroGradientFlowdDrivingStressdChin   = getDynamicColumnSizeMatrixMap< floatType, 3 * tot_dim >( d2MicroGradientFlowdDrivingStressdChin->data( ), ( num_configs - 1 ) * sot_dim );

            auto map_dPlasticMacroLdMacroStress               = dPlasticMacroLdMacroStress.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticMacroLdMicroStress               = dPlasticMacroLdMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticMicroLdMicroStress               = dPlasticMicroLdMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticGradientMicroLdMicroStress       = dPlasticGradientMicroLdMicroStress.zeroMap< floatType, tot_dim, sot_dim >( );
            auto map_dPlasticMacroLdF                         = dPlasticMacroLdF.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticMacroLdFn                        = dPlasticMacroLdFn.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticMicroLdF                         = dPlasticMicroLdF.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticMicroLdFn                        = dPlasticMicroLdFn.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticGradientMicroLdF                 = dPlasticGradientMicroLdF.zeroMap< floatType, tot_dim, sot_dim >( );
            auto map_dPlasticGradientMicroLdFn                = dPlasticGradientMicroLdFn.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticMicroLdChi                       = dPlasticMicroLdChi.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dPlasticMicroLdChin                      = dPlasticMicroLdChin.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticGradientMicroLdChi               = dPlasticGradientMicroLdChi.zeroMap< floatType, tot_dim, sot_dim >( );
            auto map_dPlasticGradientMicroLdChin              = dPlasticGradientMicroLdChin.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticGradientMicroLdGradChi           = dPlasticGradientMicroLdGradChi.zeroMap< floatType, tot_dim, tot_dim >( );
            auto map_dPlasticGradientMicroLdGradChin          = dPlasticGradientMicroLdGradChin.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * tot_dim );
            auto map_dPlasticMicroLdISVs                      = dPlasticMicroLdISVs.zeroMap< floatType, sot_dim >( num_isvs );
            auto map_dPlasticGradientMicroLdISVs              = dPlasticGradientMicroLdISVs.zeroMap< floatType, tot_dim >( num_isvs );
            auto map_dPlasticGradientMicroLdHigherOrderStress = dPlasticGradientMicroLdHigherOrderStress.zeroMap< floatType, tot_dim, tot_dim >( );

            // Assemble the Jacobians
            map_dPlasticMacroLdMacroStress = ( map_dPlasticMacroLdMacroFlowDirection * map_d2MacroFlowdDrivingStressdStress ).eval( );

            map_dPlasticMacroLdMicroStress = ( map_dPlasticMacroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdStress ).eval( );

            map_dPlasticMicroLdMicroStress = ( map_dPlasticMicroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdStress ).eval( );

            map_dPlasticGradientMicroLdMicroStress = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdMicroStress ).eval( );

            map_dPlasticGradientMicroLdHigherOrderStress = ( map_dPlasticGradientMicroLdMicroGradientFlowDirection * map_d2MicroGradientFlowdDrivingStressdStress ).eval( );

            map_dPlasticMacroLdF  = ( map_dPlasticMacroLdMacroFlowDirection * map_d2MacroFlowdDrivingStressdF ).eval( );
            map_dPlasticMacroLdF += ( map_dPlasticMacroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdF ).eval( );
            map_dPlasticMacroLdF += ( map_dPlasticMacroLdPrecedingRCG * map_dRCGdF ).eval( );

            map_dPlasticMacroLdFn  = ( map_dPlasticMacroLdMacroFlowDirection * map_d2MacroFlowdDrivingStressdFn ).eval( );
            map_dPlasticMacroLdFn += ( map_dPlasticMacroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdFn ).eval( );
            map_dPlasticMacroLdFn += ( map_dPlasticMacroLdPrecedingRCG * map_dRCGdFn ).eval( );

            map_dPlasticMicroLdF  = ( map_dPlasticMicroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdF ).eval( );
            map_dPlasticMicroLdF += ( map_dPlasticMicroLdPrecedingPsi * map_dPsidF ).eval( );

            map_dPlasticMicroLdFn  = ( map_dPlasticMicroLdMicroFlowDirection * map_d2MicroFlowdDrivingStressdFn ).eval( );
            map_dPlasticMicroLdFn += ( map_dPlasticMicroLdPrecedingPsi * map_dPsidFn ).eval( );

            map_dPlasticGradientMicroLdF  = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdF ).eval( );
            map_dPlasticGradientMicroLdF += ( map_dPlasticGradientMicroLdPrecedingPsi * map_dPsidF ).eval( );
            map_dPlasticGradientMicroLdF += ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadF ).eval( );
            map_dPlasticGradientMicroLdF += ( map_dPlasticGradientMicroLdMicroGradientFlowDirection * map_d2MicroGradientFlowdDrivingStressdF ).eval( );

            map_dPlasticGradientMicroLdFn  = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdFn ).eval( );
            map_dPlasticGradientMicroLdFn += ( map_dPlasticGradientMicroLdPrecedingPsi * map_dPsidFn ).eval( );
            map_dPlasticGradientMicroLdFn += ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadFn ).eval( );
            map_dPlasticGradientMicroLdFn += ( map_dPlasticGradientMicroLdMicroGradientFlowDirection * map_d2MicroGradientFlowdDrivingStressdFn ).eval( );

            map_dPlasticMicroLdChi  = ( map_dPlasticMicroLdPrecedingPsi * map_dPsidChi ).eval( );
            map_dPlasticMicroLdChi += ( map_dPlasticMicroLdPrecedingMicroRCG * map_dMicroRCGdChi ).eval( );

            map_dPlasticMicroLdChin  = ( map_dPlasticMicroLdPrecedingPsi * map_dPsidChin ).eval( );
            map_dPlasticMicroLdChin += ( map_dPlasticMicroLdPrecedingMicroRCG * map_dMicroRCGdChin ).eval( );

            map_dPlasticGradientMicroLdChi  = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdChi ).eval( );
            map_dPlasticGradientMicroLdChi += ( map_dPlasticGradientMicroLdPrecedingPsi * map_dPsidChi ).eval( );
            map_dPlasticGradientMicroLdChi += ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadChi ).eval( );
            map_dPlasticGradientMicroLdChi += ( map_dPlasticGradientMicroLdMicroGradientFlowDirection * map_d2MicroGradientFlowdDrivingStressdChi ).eval( );

            map_dPlasticGradientMicroLdChin  = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdChin ).eval( );
            map_dPlasticGradientMicroLdChin += ( map_dPlasticGradientMicroLdPrecedingPsi * map_dPsidChin ).eval( );
            map_dPlasticGradientMicroLdChin += ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadChin ).eval( );
            map_dPlasticGradientMicroLdChin += ( map_dPlasticGradientMicroLdMicroGradientFlowDirection * map_d2MicroGradientFlowdDrivingStressdChin ).eval( );

            map_dPlasticGradientMicroLdGradChi = ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadGradChi ).eval( );

            map_dPlasticGradientMicroLdGradChin = ( map_dPlasticGradientMicroLdPrecedingGamma * map_dPrecedingGammadGradChin ).eval( );

            dPlasticMacroLdISVs.zero( sot_dim * num_isvs );

            dPlasticMicroLdISVs.zero( sot_dim * num_isvs );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dPlasticMacroLdISVs.value )[ num_isvs * i + 0 ] = dPlasticMacroLdMacroGamma[ i ];

                ( *dPlasticMacroLdISVs.value )[ num_isvs * i + 1 ] = dPlasticMacroLdMicroGamma[ i ];

                ( *dPlasticMicroLdISVs.value )[ num_isvs * i + 1 ] = dPlasticMicroLdMicroGamma[ i ];

            }

            dPlasticGradientMicroLdISVs.zero( tot_dim * num_isvs );

            map_dPlasticGradientMicroLdISVs = ( map_dPlasticGradientMicroLdPlasticMicroL * map_dPlasticMicroLdISVs ).eval( );

            for ( unsigned int i = 0; i < tot_dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    ( *dPlasticGradientMicroLdISVs.value )[ num_isvs * i + j + 2 ] = dPlasticGradientMicroLdMicroGradientGamma[ 3 * i + j ];

                }

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

            auto updatedPlasticDeformationGradient = get_setDataStorage_updatedPlasticDeformationGradient( );

            auto updatedPlasticMicroDeformation = get_setDataStorage_updatedPlasticMicroDeformation( );

            auto updatedPlasticGradientMicroDeformation = get_setDataStorage_updatedPlasticGradientMicroDeformation( );

            const secondOrderTensor previousPlasticDeformationGradient      = secondOrderTensor( hydra->get_previousConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                                 hydra->get_previousConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const secondOrderTensor previousPlasticMicroDeformation         = secondOrderTensor( hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                                 hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const thirdOrderTensor previousPlasticGradientMicroDeformation  = thirdOrderTensor( hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
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
                                          *updatedPlasticDeformationGradient.value,
                                          *updatedPlasticMicroDeformation.value,
                                          *updatedPlasticGradientMicroDeformation.value,
                                          *getIntegrationParameter( ),
                                          *getIntegrationParameter( ),
                                          *getIntegrationParameter( ) );
            )

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

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            auto updatedPlasticDeformationGradient = get_setDataStorage_updatedPlasticDeformationGradient( );

            auto updatedPlasticMicroDeformation = get_setDataStorage_updatedPlasticMicroDeformation( );

            auto updatedPlasticGradientMicroDeformation = get_setDataStorage_updatedPlasticGradientMicroDeformation( );

            const secondOrderTensor previousPlasticDeformationGradient      = secondOrderTensor( hydra->get_previousConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                                 hydra->get_previousConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const secondOrderTensor previousPlasticMicroDeformation         = secondOrderTensor( hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                                 hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const thirdOrderTensor  previousPlasticGradientMicroDeformation = thirdOrderTensor( hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
                                                                                                hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim * ( plasticConfigurationIndex + 1 ) );

            fourthOrderTensor dPlasticFdPlasticMacroL;

            fourthOrderTensor dPlasticMicroDeformationdPlasticMicroL;

            fifthOrderTensor  dPlasticGradientMicroDeformationdPlasticMacroL;

            fifthOrderTensor  dPlasticGradientMicroDeformationdPlasticMicroL;

            sixthOrderTensor  dPlasticGradientMicroDeformationdPlasticGradientMicroL;

            if ( addPreviousGradients ){

                fourthOrderTensor dPlasticFdPreviousPlasticF;
                fourthOrderTensor dPlasticFdPreviousPlasticMacroL;
                fourthOrderTensor dPlasticMicroDeformationdPreviousPlasticMicroDeformation;
                fourthOrderTensor dPlasticMicroDeformationdPreviousPlasticMicroL;
                fifthOrderTensor  dPlasticGradientMicroDeformationdPreviousPlasticMicroDeformation;
                sixthOrderTensor  dPlasticGradientMicroDeformationdPreviousPlasticMicroGradient;
                fifthOrderTensor  dPlasticGradientMicroDeformationdPreviousPlasticMacroL;
                fifthOrderTensor  dPlasticGradientMicroDeformationdPreviousPlasticMicroL;
                sixthOrderTensor  dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL;

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
                                              *updatedPlasticDeformationGradient.value,
                                              *updatedPlasticMicroDeformation.value,
                                              *updatedPlasticGradientMicroDeformation.value,
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

                auto map_dPlasticFdPreviousPlasticMacroL                                  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticFdPreviousPlasticMacroL.data( ) );
                auto map_dPlasticMicroDeformationdPreviousPlasticMicroL                   = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMicroDeformationdPreviousPlasticMicroL.data( ) );
                auto map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL           = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPlasticGradientMicroDeformationdPreviousPlasticMacroL.data( ) );
                auto map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL           = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPlasticGradientMicroDeformationdPreviousPlasticMicroL.data( ) );
                auto map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL   = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL.data( ) );

                auto map_previousdPlasticMacroVelocityGradientdMacroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMacroVelocityGradientdMacroStress( )->data( ) );
                auto map_previousdPlasticMacroVelocityGradientdMicroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMacroVelocityGradientdMicroStress( )->data( ) );
                auto map_previousdPlasticMacroVelocityGradientdF                          = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMacroVelocityGradientdF( )->data( ) );
                auto map_previousdPlasticMacroVelocityGradientdFn                         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_previousdPlasticMacroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticMacroVelocityGradientdStateVariables             = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_previousdPlasticMacroVelocityGradientdStateVariables( )->data( ), num_isvs );
                auto map_previousdPlasticMicroVelocityGradientdMicroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMicroVelocityGradientdMicroStress( )->data( ) );
                auto map_previousdPlasticMicroVelocityGradientdF                          = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMicroVelocityGradientdF( )->data( ) );
                auto map_previousdPlasticMicroVelocityGradientdFn                         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_previousdPlasticMicroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticMicroVelocityGradientdChi                        = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_previousdPlasticMicroVelocityGradientdChi( )->data( ) );
                auto map_previousdPlasticMicroVelocityGradientdChin                       = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_previousdPlasticMicroVelocityGradientdChin( )->data( ), sot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticMicroVelocityGradientdStateVariables             = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_previousdPlasticMicroVelocityGradientdStateVariables( )->data( ), num_isvs );
                auto map_previousdPlasticGradientMicroVelocityGradientdMicroStress        = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_previousdPlasticGradientMicroVelocityGradientdMicroStress( )->data( ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress  = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( )->data( ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdF                  = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_previousdPlasticGradientMicroVelocityGradientdF( )->data( ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdFn                 = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdChi                = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_previousdPlasticGradientMicroVelocityGradientdChi( )->data( ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdChin               = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdChin( )->data( ), sot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdGradChi            = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdGradChi( )->data( ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdGradChin           = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdGradChin( )->data( ), tot_dim * ( num_configs - 1 ) );
                auto map_previousdPlasticGradientMicroVelocityGradientdStateVariables     = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_previousdPlasticGradientMicroVelocityGradientdStateVariables( )->data( ), num_isvs );

                auto dUpdatedPlasticDeformationGradientdPreviousMacroStress = get_setDataStorage_dUpdatedPlasticDeformationGradientdPreviousMacroStress( );
                auto map_dUpdatedPlasticDeformationGradientdPreviousMacroStress = dUpdatedPlasticDeformationGradientdPreviousMacroStress.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dUpdatedPlasticDeformationGradientdPreviousMicroStress = get_setDataStorage_dUpdatedPlasticDeformationGradientdPreviousMicroStress( );
                auto map_dUpdatedPlasticDeformationGradientdPreviousMicroStress = dUpdatedPlasticDeformationGradientdPreviousMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dUpdatedPlasticDeformationGradientdPreviousF = get_setDataStorage_dUpdatedPlasticDeformationGradientdPreviousF( );
                auto map_dUpdatedPlasticDeformationGradientdPreviousF = dUpdatedPlasticDeformationGradientdPreviousF.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dUpdatedPlasticFdPreviousFn = get_setDataStorage_dUpdatedPlasticDeformationGradientdPreviousFn( );
                auto map_dUpdatedPlasticFdPreviousFn = dUpdatedPlasticFdPreviousFn.zeroMap< floatType, sot_dim >( sot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticDeformationGradientdPreviousStateVariables = get_setDataStorage_dUpdatedPlasticDeformationGradientdPreviousStateVariables( );
                auto map_dUpdatedPlasticDeformationGradientdPreviousStateVariables = dUpdatedPlasticDeformationGradientdPreviousStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

                auto dUpdatedPlasticMicroDeformationdPreviousMicroStress = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousMicroStress( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousMicroStress = dUpdatedPlasticMicroDeformationdPreviousMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dUpdatedPlasticMicroDeformationdPreviousF = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousF( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousF = dUpdatedPlasticMicroDeformationdPreviousF.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dUpdatedPlasticMicroDeformationdPreviousFn = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousFn( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousFn = dUpdatedPlasticMicroDeformationdPreviousFn.zeroMap< floatType, sot_dim >( sot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticMicroDeformationdPreviousChi = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousChi( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousChi = dUpdatedPlasticMicroDeformationdPreviousChi.zeroMap< floatType, sot_dim >( sot_dim );

                auto dUpdatedPlasticMicroDeformationdPreviousChin = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousChin( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousChin = dUpdatedPlasticMicroDeformationdPreviousChin.zeroMap< floatType, sot_dim >( sot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticMicroDeformationdPreviousStateVariables = get_setDataStorage_dUpdatedPlasticMicroDeformationdPreviousStateVariables( );
                auto map_dUpdatedPlasticMicroDeformationdPreviousStateVariables = dUpdatedPlasticMicroDeformationdPreviousStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress = dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress.zeroMap< floatType, tot_dim, sot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress = dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress.zeroMap< floatType, tot_dim, sot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress = dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress.zeroMap< floatType, tot_dim, tot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousF = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousF( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousF = dUpdatedPlasticGradientMicroDeformationdPreviousF.zeroMap< floatType, tot_dim, sot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousFn = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousFn( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousFn = dUpdatedPlasticGradientMicroDeformationdPreviousFn.zeroMap< floatType, tot_dim >( sot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousChi = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousChi( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousChi = dUpdatedPlasticGradientMicroDeformationdPreviousChi.zeroMap< floatType, tot_dim, sot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousChin = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousChin( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousChin = dUpdatedPlasticGradientMicroDeformationdPreviousChin.zeroMap< floatType, tot_dim >( sot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousGradChi = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi = dUpdatedPlasticGradientMicroDeformationdPreviousGradChi.zeroMap< floatType, tot_dim, tot_dim >( );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousGradChin = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin = dUpdatedPlasticGradientMicroDeformationdPreviousGradChin.zeroMap< floatType, tot_dim >( tot_dim * ( num_configs - 1 ) );

                auto dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables( );
                auto map_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables = dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables.zeroMap< floatType, tot_dim >( num_isvs );

                map_dUpdatedPlasticDeformationGradientdPreviousMacroStress = ( map_dPlasticFdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdMacroStress ).eval( );

                map_dUpdatedPlasticDeformationGradientdPreviousMicroStress = ( map_dPlasticFdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdMicroStress );

                map_dUpdatedPlasticDeformationGradientdPreviousF = ( map_dPlasticFdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdF ).eval( );

                map_dUpdatedPlasticFdPreviousFn = ( map_dPlasticFdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdFn ).eval( );

                unsigned int offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dUpdatedPlasticFdPreviousFn.value )[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticFdPreviousPlasticF[ sot_dim * i + j ];

                    }

                }

                map_dUpdatedPlasticDeformationGradientdPreviousStateVariables = ( map_dPlasticFdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdStateVariables ).eval( );

                map_dUpdatedPlasticMicroDeformationdPreviousMicroStress = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdMicroStress ).eval( );

                map_dUpdatedPlasticMicroDeformationdPreviousF = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdF ).eval( );

                map_dUpdatedPlasticMicroDeformationdPreviousFn = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdFn ).eval( );

                map_dUpdatedPlasticMicroDeformationdPreviousChi = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdChi ).eval( );

                map_dUpdatedPlasticMicroDeformationdPreviousChin = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdChin ).eval( );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dUpdatedPlasticMicroDeformationdPreviousChin.value )[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticMicroDeformationdPreviousPlasticMicroDeformation[ sot_dim * i + j ];

                    }

                }

                map_dUpdatedPlasticMicroDeformationdPreviousStateVariables = ( map_dPlasticMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdStateVariables ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdMacroStress ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdMicroStress ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress += ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdMicroStress ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdMicroStress ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress = ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousF  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdF );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousF += ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdF );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousF += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdF );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousFn  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdFn ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousFn += ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdFn ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousFn += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdFn ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousChi  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdChi ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousChi += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdChi ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousChin  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdChin ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousChin += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdChin ).eval( );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * sot_dim;

                for ( unsigned int i = 0; i < tot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dUpdatedPlasticGradientMicroDeformationdPreviousChin.value )[ ( num_configs - 1 ) * sot_dim * i + j + offset ] += dPlasticGradientMicroDeformationdPreviousPlasticMicroDeformation[ sot_dim * i + j ];

                    }

                }

                map_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi = ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdGradChi ).eval( );

                map_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin = ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdGradChin ).eval( );

                offset = ( ( *getPlasticConfigurationIndex( ) ) - 1 ) * tot_dim;

                for ( unsigned int i = 0; i < tot_dim; i++ ){

                    for ( unsigned int j = 0; j < tot_dim; j++ ){

                        ( *dUpdatedPlasticGradientMicroDeformationdPreviousGradChin.value )[ ( num_configs - 1 ) * tot_dim * i + j + offset ] += dPlasticGradientMicroDeformationdPreviousPlasticMicroGradient[ tot_dim * i  + j ];

                    }

                }

                map_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables  = ( map_dPlasticGradientMicroDeformationdPreviousPlasticMacroL * map_previousdPlasticMacroVelocityGradientdStateVariables ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables += ( map_dPlasticGradientMicroDeformationdPreviousPlasticMicroL * map_previousdPlasticMicroVelocityGradientdStateVariables ).eval( );
                map_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables += ( map_dPlasticGradientMicroDeformationdPreviousPlasticGradientMicroL * map_previousdPlasticGradientMicroVelocityGradientdStateVariables ).eval( );

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
                                              *updatedPlasticDeformationGradient.value,
                                              *updatedPlasticMicroDeformation.value,
                                              *updatedPlasticGradientMicroDeformation.value,
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

            auto map_dPlasticFdPlasticMacroL                                  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticFdPlasticMacroL.data( ) );
            auto map_dPlasticMicroDeformationdPlasticMicroL                   = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dPlasticMicroDeformationdPlasticMicroL.data( ) );
            auto map_dPlasticGradientMicroDeformationdPlasticMacroL           = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPlasticGradientMicroDeformationdPlasticMacroL.data( ) );
            auto map_dPlasticGradientMicroDeformationdPlasticMicroL           = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( dPlasticGradientMicroDeformationdPlasticMicroL.data( ) );
            auto map_dPlasticGradientMicroDeformationdPlasticGradientMicroL   = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( dPlasticGradientMicroDeformationdPlasticGradientMicroL.data( ) );

            auto map_dPlasticMacroVelocityGradientdMacroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMacroVelocityGradientdMacroStress( )->data( ) );
            auto map_dPlasticMacroVelocityGradientdMicroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMacroVelocityGradientdMicroStress( )->data( ) );
            auto map_dPlasticMacroVelocityGradientdF                          = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMacroVelocityGradientdF( )->data( ) );
            auto map_dPlasticMacroVelocityGradientdFn                         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPlasticMacroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
            auto map_dPlasticMacroVelocityGradientdStateVariables             = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPlasticMacroVelocityGradientdStateVariables( )->data( ), num_isvs );
            auto map_dPlasticMicroVelocityGradientdMicroStress                = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMicroVelocityGradientdMicroStress( )->data( ) );
            auto map_dPlasticMicroVelocityGradientdF                          = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMicroVelocityGradientdF( )->data( ) );
            auto map_dPlasticMicroVelocityGradientdFn                         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPlasticMicroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
            auto map_dPlasticMicroVelocityGradientdChi                        = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPlasticMicroVelocityGradientdChi( )->data( ) );
            auto map_dPlasticMicroVelocityGradientdChin                       = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPlasticMicroVelocityGradientdChin( )->data( ), sot_dim * ( num_configs - 1 ) );
            auto map_dPlasticMicroVelocityGradientdStateVariables             = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPlasticMicroVelocityGradientdStateVariables( )->data( ), num_isvs );
            auto map_dPlasticGradientMicroVelocityGradientdMicroStress        = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_dPlasticGradientMicroVelocityGradientdMicroStress( )->data( ) );
            auto map_dPlasticGradientMicroVelocityGradientdHigherOrderStress  = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( get_dPlasticGradientMicroVelocityGradientdHigherOrderStress( )->data( ) );
            auto map_dPlasticGradientMicroVelocityGradientdF                  = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_dPlasticGradientMicroVelocityGradientdF( )->data( ) );
            auto map_dPlasticGradientMicroVelocityGradientdFn                 = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_dPlasticGradientMicroVelocityGradientdFn( )->data( ), sot_dim * ( num_configs - 1 ) );
            auto map_dPlasticGradientMicroVelocityGradientdChi                = getFixedSizeMatrixMap< floatType, tot_dim, sot_dim >( get_dPlasticGradientMicroVelocityGradientdChi( )->data( ) );
            auto map_dPlasticGradientMicroVelocityGradientdChin               = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_dPlasticGradientMicroVelocityGradientdChin( )->data( ), sot_dim * ( num_configs - 1 ) );
            auto map_dPlasticGradientMicroVelocityGradientdGradChi            = getFixedSizeMatrixMap< floatType, tot_dim, tot_dim >( get_dPlasticGradientMicroVelocityGradientdGradChi( )->data( ) );
            auto map_dPlasticGradientMicroVelocityGradientdGradChin           = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_dPlasticGradientMicroVelocityGradientdGradChin( )->data( ), tot_dim * ( num_configs - 1 ) );
            auto map_dPlasticGradientMicroVelocityGradientdStateVariables     = getDynamicColumnSizeMatrixMap< floatType, tot_dim >( get_dPlasticGradientMicroVelocityGradientdStateVariables( )->data( ), num_isvs );

            auto dUpdatedPlasticDeformationGradientdMacroStress = get_setDataStorage_dUpdatedPlasticDeformationGradientdMacroStress( );
            auto map_dUpdatedPlasticDeformationGradientdMacroStress = dUpdatedPlasticDeformationGradientdMacroStress.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticDeformationGradientdMicroStress = get_setDataStorage_dUpdatedPlasticDeformationGradientdMicroStress( );
            auto map_dUpdatedPlasticDeformationGradientdMicroStress = dUpdatedPlasticDeformationGradientdMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticDeformationGradientdF = get_setDataStorage_dUpdatedPlasticDeformationGradientdF( );
            auto map_dUpdatedPlasticDeformationGradientdF = dUpdatedPlasticDeformationGradientdF.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticDeformationGradientdFn = get_setDataStorage_dUpdatedPlasticDeformationGradientdFn( );
            auto map_dUpdatedPlasticDeformationGradientdFn = dUpdatedPlasticDeformationGradientdFn.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dUpdatedPlasticDeformationGradientdStateVariables = get_setDataStorage_dUpdatedPlasticDeformationGradientdStateVariables( );
            auto map_dUpdatedPlasticDeformationGradientdStateVariables = dUpdatedPlasticDeformationGradientdStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

            auto dUpdatedPlasticMicroDeformationdMicroStress = get_setDataStorage_dUpdatedPlasticMicroDeformationdMicroStress( );
            auto map_dUpdatedPlasticMicroDeformationdMicroStress = dUpdatedPlasticMicroDeformationdMicroStress.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticMicroDeformationdF = get_setDataStorage_dUpdatedPlasticMicroDeformationdF( );
            auto map_dUpdatedPlasticMicroDeformationdF = dUpdatedPlasticMicroDeformationdF.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticMicroDeformationdFn = get_setDataStorage_dUpdatedPlasticMicroDeformationdFn( );
            auto map_dUpdatedPlasticMicroDeformationdFn = dUpdatedPlasticMicroDeformationdFn.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dUpdatedPlasticMicroDeformationdChi = get_setDataStorage_dUpdatedPlasticMicroDeformationdChi( );
            auto map_dUpdatedPlasticMicroDeformationdChi = dUpdatedPlasticMicroDeformationdChi.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dUpdatedPlasticMicroDeformationdChin = get_setDataStorage_dUpdatedPlasticMicroDeformationdChin( );
            auto map_dUpdatedPlasticMicroDeformationdChin = dUpdatedPlasticMicroDeformationdChin.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dUpdatedPlasticMicroDeformationdStateVariables = get_setDataStorage_dUpdatedPlasticMicroDeformationdStateVariables( );
            auto map_dUpdatedPlasticMicroDeformationdStateVariables = dUpdatedPlasticMicroDeformationdStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

            auto dUpdatedPlasticGradientMicroDeformationdMacroStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdMacroStress( );
            auto map_dUpdatedPlasticGradientMicroDeformationdMacroStress = dUpdatedPlasticGradientMicroDeformationdMacroStress.zeroMap< floatType, tot_dim, sot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdMicroStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdMicroStress( );
            auto map_dUpdatedPlasticGradientMicroDeformationdMicroStress = dUpdatedPlasticGradientMicroDeformationdMicroStress.zeroMap< floatType, tot_dim, sot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdHigherOrderStress = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress( );
            auto map_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress = dUpdatedPlasticGradientMicroDeformationdHigherOrderStress.zeroMap< floatType, tot_dim, tot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdF = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdF( );
            auto map_dUpdatedPlasticGradientMicroDeformationdF = dUpdatedPlasticGradientMicroDeformationdF.zeroMap< floatType, tot_dim, sot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdFn = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdFn( );
            auto map_dUpdatedPlasticGradientMicroDeformationdFn = dUpdatedPlasticGradientMicroDeformationdFn.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dUpdatedPlasticGradientMicroDeformationdChi = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdChi( );
            auto map_dUpdatedPlasticGradientMicroDeformationdChi = dUpdatedPlasticGradientMicroDeformationdChi.zeroMap< floatType, tot_dim, sot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdChin = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdChin( );
            auto map_dUpdatedPlasticGradientMicroDeformationdChin = dUpdatedPlasticGradientMicroDeformationdChin.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dUpdatedPlasticGradientMicroDeformationdGradChi = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdGradChi( );
            auto map_dUpdatedPlasticGradientMicroDeformationdGradChi = dUpdatedPlasticGradientMicroDeformationdGradChi.zeroMap< floatType, tot_dim, tot_dim >( );

            auto dUpdatedPlasticGradientMicroDeformationdGradChin = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdGradChin( );
            auto map_dUpdatedPlasticGradientMicroDeformationdGradChin = dUpdatedPlasticGradientMicroDeformationdGradChin.zeroMap< floatType, tot_dim >( ( num_configs - 1 ) * tot_dim );

            auto dUpdatedPlasticGradientMicroDeformationdStateVariables = get_setDataStorage_dUpdatedPlasticGradientMicroDeformationdStateVariables( );
            auto map_dUpdatedPlasticGradientMicroDeformationdStateVariables = dUpdatedPlasticGradientMicroDeformationdStateVariables.zeroMap< floatType, tot_dim >( num_isvs );

            map_dUpdatedPlasticDeformationGradientdMacroStress = ( map_dPlasticFdPlasticMacroL * map_dPlasticMacroVelocityGradientdMacroStress ).eval( );

            map_dUpdatedPlasticDeformationGradientdMicroStress = ( map_dPlasticFdPlasticMacroL * map_dPlasticMacroVelocityGradientdMicroStress ).eval( );

            map_dUpdatedPlasticDeformationGradientdF = ( map_dPlasticFdPlasticMacroL * map_dPlasticMacroVelocityGradientdF ).eval( );

            map_dUpdatedPlasticDeformationGradientdFn = ( map_dPlasticFdPlasticMacroL * map_dPlasticMacroVelocityGradientdFn ).eval( );

            map_dUpdatedPlasticDeformationGradientdStateVariables = ( map_dPlasticFdPlasticMacroL * map_dPlasticMacroVelocityGradientdStateVariables ).eval( );

            map_dUpdatedPlasticMicroDeformationdMicroStress = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdMicroStress ).eval( );

            map_dUpdatedPlasticMicroDeformationdF = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdF ).eval( );

            map_dUpdatedPlasticMicroDeformationdFn = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdFn ).eval( );        

            map_dUpdatedPlasticMicroDeformationdChi = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdChi ).eval( );

            map_dUpdatedPlasticMicroDeformationdChin = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdChin ).eval( );

            map_dUpdatedPlasticMicroDeformationdStateVariables = ( map_dPlasticMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdStateVariables ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdMacroStress = ( map_dPlasticGradientMicroDeformationdPlasticMacroL * map_dPlasticMacroVelocityGradientdMacroStress ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdMicroStress  = ( map_dPlasticGradientMicroDeformationdPlasticMacroL * map_dPlasticMacroVelocityGradientdMicroStress ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdMicroStress += ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdMicroStress ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdMicroStress += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdMicroStress ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress = ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdHigherOrderStress ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdF  = ( map_dPlasticGradientMicroDeformationdPlasticMacroL * map_dPlasticMacroVelocityGradientdF ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdF += ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdF ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdF += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdF ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdFn  = ( map_dPlasticGradientMicroDeformationdPlasticMacroL * map_dPlasticMacroVelocityGradientdFn ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdFn += ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdFn ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdFn += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdFn ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdChi  = ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdChi ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdChi += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdChi ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdChin  = ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdChin ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdChin += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdChin ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdGradChi = ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdGradChi ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdGradChin = ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdGradChin ).eval( );

            map_dUpdatedPlasticGradientMicroDeformationdStateVariables  = ( map_dPlasticGradientMicroDeformationdPlasticMacroL * map_dPlasticMacroVelocityGradientdStateVariables ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdStateVariables += ( map_dPlasticGradientMicroDeformationdPlasticMicroL * map_dPlasticMicroVelocityGradientdStateVariables ).eval( );
            map_dUpdatedPlasticGradientMicroDeformationdStateVariables += ( map_dPlasticGradientMicroDeformationdPlasticGradientMicroL * map_dPlasticGradientMicroVelocityGradientdStateVariables ).eval( );

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

            const dimVector *microGradientYield = get_microGradientYield( );

            auto residual = get_setDataStorage_stateVariableResiduals( );
            residual.zero( get_plasticStateVariables( )->size( ) );

            floatType macroMac  = tardigradeConstitutiveTools::mac( *macroYield );
            tardigradeConstitutiveTools::mac( -( *macroYield ) );

            floatType microMac  = tardigradeConstitutiveTools::mac( *microYield );
            tardigradeConstitutiveTools::mac( -( *microYield ) );

            floatVector microGradientMac( microGradientYield->size( ), 0 );

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
                    macNegMicroGradientGamma[ index ] = tardigradeConstitutiveTools::mac( -( *plasticMultipliers )[ index + 2 ] );
                }

            }

            // Set the terms associated with the yield surface
            ( *residual.value )[ 0 ] = macroMac + ( *plasticMultipliers )[ 0 ] * ( *macroYield ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMacroGamma;

            ( *residual.value )[ 1 ] = microMac + ( *plasticMultipliers )[ 1 ] * ( *microYield ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMicroGamma;

            for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){

                unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );

                ( *residual.value )[ index + 2 ]
                    = microGradientMac[ index ] + ( *plasticMultipliers )[ index + 2 ] * ( *y ) + ( *getPlasticMultiplierBarrierModulus( ) ) * macNegMicroGradientGamma[ index ];

            }

            // Set the terms associated with the strain-like ISV evolution
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *residual.value )[ numPlasticMultipliers + i ] = ( *updatedPlasticStrainLikeISVs )[ i ] - ( *plasticStrainLikeISVs )[ i ];

            }

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

            const unsigned int numUnknowns = hydra->getNumUnknowns( );

            const unsigned int numISVs = get_plasticStateVariables( )->size( );

            const secondOrderTensor *dMacroYielddStress                 = get_dMacroYielddStress( );

            const floatVector *dMacroYielddFn                     = get_dMacroYielddFn( );

            const floatVector *dMacroYielddStateVariables         = get_dMacroYielddStateVariables( );

            const secondOrderTensor *dMicroYielddStress                 = get_dMicroYielddStress( );

            const floatVector *dMicroYielddFn                     = get_dMicroYielddFn( );

            const floatVector *dMicroYielddStateVariables         = get_dMicroYielddStateVariables( );

            const floatVector *dMicroGradientYielddStress         = get_dMicroGradientYielddStress( );

            const floatVector *dMicroGradientYielddFn             = get_dMicroGradientYielddFn( );

            const floatVector *dMicroGradientYielddChin           = get_dMicroGradientYielddChin( );

            const floatVector *dMicroGradientYielddStateVariables = get_dMicroGradientYielddStateVariables( );

            const floatVector *dUpdatedPlasticStrainLikeISVsdStateVariables = get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto jacobian = get_setDataStorage_stateVariableJacobians( );
            jacobian.zero( numISVs * numUnknowns );

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

                ( *jacobian.value )[ numUnknowns * 0 + j ] = ( dMacroMacdx  + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddStress )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddStress )[ j ];

            }

            offset = 2 * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddStress )[ numThirdOrderTensor * i + j ];

                }

            }

            // Sub-Deformation gradient jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor;
            for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                ( *jacobian.value )[ numUnknowns * 0 + j + offset ] = ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddFn )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddFn )[ j ];

            } 

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }

            // Sub-Micro deformation jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor + ( numConfigurations - 1 ) * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }


            // State Variable Jacobians
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );

            ( *jacobian.value )[ numUnknowns * 0 + offset + 0 ] += *macroYield - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMacroGammadGamma;

            ( *jacobian.value )[ numUnknowns * 1 + offset + 1 ] += *microYield - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMicroGammadGamma;

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + 2 ) + offset + i + 2 ] += ( *microGradientYield )[ i ] - ( *getPlasticMultiplierBarrierModulus( ) ) * dMacNegMicroGradientGammadGamma[ i ];

            }

            for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                ( *jacobian.value )[ numUnknowns * 0 + j + offset ] += ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddStateVariables )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] += ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddStateVariables )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] += ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddStateVariables )[ numISVs * i + j ];

                }

            }
            
            unsigned int row0 = numPlasticMultipliers;
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset + numPlasticMultipliers ] -= 1;

                for ( auto j = getStateVariableIndices( )->begin( ); j != getStateVariableIndices( )->end( ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numISVs * i + ( unsigned int )( j - getStateVariableIndices( )->begin( ) ) ];

                }

            }

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

            const secondOrderTensor *dMacroYielddF                      = get_dMacroYielddF( );

            const secondOrderTensor *dMicroYielddF                      = get_dMicroYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddF              = get_dMicroGradientYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddChi            = get_dMicroGradientYielddChi( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto dRdD = get_setDataStorage_dStateVariableResidualsdD( );
            dRdD.zero( get_plasticStateVariables( )->size( ) * numConfigurationUnknowns );

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

                ( *dRdD.value )[ numConfigurationUnknowns * 0 + j + offset ] = ( dMacroMacdx + ( *plasticMultipliers )[ 0 ] ) * ( *dMacroYielddF )[ j ];

                ( *dRdD.value )[ numConfigurationUnknowns * 1 + j + offset ] = ( dMicroMacdx + ( *plasticMultipliers )[ 1 ] ) * ( *dMicroYielddF )[ j ];

            } 

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddF )[ numSecondOrderTensor * i + j ];

                }

            }

            // Micro deformation jacobians
            offset = numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( dMicroGradientMacdx[ i ] + ( *plasticMultipliers )[ i + 2 ] ) * ( *dMicroGradientYielddChi )[ numSecondOrderTensor * i + j ];

                }

            }

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

            auto dRdPreviousISVs = get_setDataStorage_dStateVariableResidualsdPreviousISVs( );
            dRdPreviousISVs.zero( numPlasticISVs * numISVs );

            // Stress Jacobians
            unsigned int row0 = numPlasticMultipliers;
            unsigned int offset = ( numConfigurations - 1 ) * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            std::vector< unsigned int > stateVariableIndices = *getStateVariableIndices( );

            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *dRdPreviousISVs.value )[ numISVs * ( i + row0 ) + stateVariableIndices[ i ] + offset + numPlasticMultipliers ] += 1;

                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    ( *dRdPreviousISVs.value )[ numISVs * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numPlasticISVs * i + ( unsigned int )( j - stateVariableIndices.begin( ) ) ];

                }

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the residual equation
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );
 
            const unsigned int tot_dim = hydra->getTOTDimension( ); 

            const secondOrderTensor *updatedPlasticDeformationGradient;

            const secondOrderTensor *updatedPlasticMicroDeformation;

            const thirdOrderTensor  *updatedPlasticGradientMicroDeformation;

            const floatVector *stateVariableResiduals;

            // Get the trial plastic deformation measures
            unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            const secondOrderTensor plasticDeformationGradient      = secondOrderTensor( hydra->get_configurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                         hydra->get_configurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const secondOrderTensor plasticMicroDeformation         = secondOrderTensor( hydra->get_microConfigurations( )->begin( ) + sot_dim * plasticConfigurationIndex,
                                                                                         hydra->get_microConfigurations( )->begin( ) + sot_dim * ( plasticConfigurationIndex + 1 ) );

            const thirdOrderTensor plasticGradientMicroDeformation  = thirdOrderTensor( hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim * plasticConfigurationIndex,
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

            auto residual = get_setDataStorage_residual( );

            *residual.value = tardigradeVectorTools::appendVectors( { *updatedPlasticDeformationGradient      - plasticDeformationGradient,
                                                                      *updatedPlasticMicroDeformation         - plasticMicroDeformation,
                                                                      *updatedPlasticGradientMicroDeformation - plasticGradientMicroDeformation,
                                                                      *stateVariableResiduals } );

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian
             */

            const unsigned int numEquations = *getNumEquations( );

            const unsigned int numUnknowns  = hydra->getNumUnknowns( );

            const unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            const std::vector< unsigned int > stateVariableIndices = *getStateVariableIndices( );

            const unsigned int numISVs = stateVariableIndices.size( );

            const fourthOrderTensor *dUpdatedPlasticDeformationGradientdMacroStress;

            const fourthOrderTensor *dUpdatedPlasticDeformationGradientdMicroStress;

            const floatVector *dUpdatedPlasticDeformationGradientdFn;

            const floatVector *dUpdatedPlasticDeformationGradientdStateVariables;

            const fourthOrderTensor *dUpdatedPlasticMicroDeformationdMicroStress;

            const floatVector *dUpdatedPlasticMicroDeformationdFn;

            const floatVector *dUpdatedPlasticMicroDeformationdChin;

            const floatVector *dUpdatedPlasticMicroDeformationdStateVariables;

            const fifthOrderTensor *dUpdatedPlasticGradientMicroDeformationdMacroStress;

            const fifthOrderTensor *dUpdatedPlasticGradientMicroDeformationdMicroStress;

            const sixthOrderTensor *dUpdatedPlasticGradientMicroDeformationdHigherOrderStress;

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

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( numEquations * numUnknowns );

            // Set the Jacobians of the second order plastic deformation measures
            for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                // Jacobians with respect to the trial stresses
                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i                        ) + j                        ] += ( *dUpdatedPlasticDeformationGradientdMacroStress )[ numSecondOrderTensor * i + j ];

                    ( *jacobian.value )[ numUnknowns * ( i                        ) + j + numSecondOrderTensor ] += ( *dUpdatedPlasticDeformationGradientdMicroStress )[ numSecondOrderTensor * i + j ];

                    ( *jacobian.value )[ numUnknowns * ( i + numSecondOrderTensor ) + j + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdMicroStress )[ numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-deformation gradients
                ( *jacobian.value )[ numUnknowns * ( i                        ) + i + 2 * numSecondOrderTensor + numThirdOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i                        ) + j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticDeformationGradientdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                    ( *jacobian.value )[ numUnknowns * ( i + numSecondOrderTensor ) + j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-micro deformation

                ( *jacobian.value )[ numUnknowns * ( i + numSecondOrderTensor ) + i + numConfigurationUnknowns + numSecondOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + numSecondOrderTensor ) + j + numConfigurationUnknowns + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the state variables
                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    unsigned int col = ( unsigned int )( j - stateVariableIndices.begin( ) );

                    ( *jacobian.value )[ numUnknowns * ( i                        ) + *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticDeformationGradientdStateVariables )[ numISVs * i + col ];

                    ( *jacobian.value )[ numUnknowns * ( i + numSecondOrderTensor ) + *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticMicroDeformationdStateVariables )[ numISVs * i + col ];

                }

            }

            // Set the Jacobians of the third order plastic deformation measures
            for ( unsigned int i = 0; i < numThirdOrderTensor; i++ ){

                // Set the Jacobians with respect to the trial stresses
                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j                        ] = ( *dUpdatedPlasticGradientMicroDeformationdMacroStress )[ numSecondOrderTensor * i + j ];

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j + numSecondOrderTensor ] = ( *dUpdatedPlasticGradientMicroDeformationdMicroStress )[ numSecondOrderTensor * i + j ];

                }

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j + 2 * numSecondOrderTensor ] = ( *dUpdatedPlasticGradientMicroDeformationdHigherOrderStress )[ numThirdOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j + 2 * numSecondOrderTensor + numThirdOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-micro deformations
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j + numConfigurationUnknowns + numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

                // Set the jacobians with respect to the local spatial gradients of the sub-micro deformation
                ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + i + numConfigurationUnknowns + 2 * numSecondOrderTensor ] -= 1;

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numThirdOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + j + numConfigurationUnknowns + 2 * numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdGradChin )[ ( numConfigurations - 1 ) * numThirdOrderTensor * i + j ];

                }

                // Jacobians with respect to the state variables
                for ( auto j = stateVariableIndices.begin( ); j != stateVariableIndices.end( ); j++ ){

                    unsigned int col = ( unsigned int )( j - stateVariableIndices.begin( ) );

                    ( *jacobian.value )[ numUnknowns * ( i + 2 * numSecondOrderTensor ) + *j + numConfigurations * numConfigurationUnknowns ] += ( *dUpdatedPlasticGradientMicroDeformationdStateVariables )[ numISVs * i + col ];

                }

            }

            // Set the Jacobians of the state variables
            std::copy( stateVariableJacobians->begin( ), stateVariableJacobians->end( ), std::begin( *jacobian.value ) + numUnknowns * numConfigurationUnknowns );

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

            const fourthOrderTensor *dUpdatedPlasticDeformationGradientdF;

            const fourthOrderTensor *dUpdatedPlasticMicroDeformationdF;

            const fourthOrderTensor *dUpdatedPlasticMicroDeformationdChi;

            const fifthOrderTensor *dUpdatedPlasticGradientMicroDeformationdF;

            const fifthOrderTensor *dUpdatedPlasticGradientMicroDeformationdChi;

            const sixthOrderTensor *dUpdatedPlasticGradientMicroDeformationdGradChi;

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

            auto dRdD = get_setDataStorage_dRdD( );
            dRdD.zero( numEquations * numConfigurationUnknowns );

            // Set the Jacobians of the second order plastic deformation measures
            for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                // Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i                        ) + j ] += ( *dUpdatedPlasticDeformationGradientdF )[ numSecondOrderTensor * i + j ];

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + numSecondOrderTensor ) + j ] += ( *dUpdatedPlasticMicroDeformationdF )[ numSecondOrderTensor * i + j ];

                }

                // Jacobians with respect to the trial sub-micro deformation
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + numSecondOrderTensor ) + j + numSecondOrderTensor ] += ( *dUpdatedPlasticMicroDeformationdChi )[ numSecondOrderTensor * i + j ];
 
                }

            }

            // Set the Jacobians of the third order plastic deformation measures
            for ( unsigned int i = 0; i < numThirdOrderTensor; i++ ){

                // Set the Jacobians with respect to the trial sub-deformation gradients
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 * numSecondOrderTensor ) + j ] += ( *dUpdatedPlasticGradientMicroDeformationdF )[ numSecondOrderTensor * i + j ];

                }

                // Set the Jacobians with respect to the trial sub-micro deformations
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 * numSecondOrderTensor ) + j + numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdChi )[ numSecondOrderTensor * i + j ];

                }

                // Set the jacobians with respect to the local spatial gradients of the sub-micro deformation
                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numThirdOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 * numSecondOrderTensor ) + j + 2 * numSecondOrderTensor ] += ( *dUpdatedPlasticGradientMicroDeformationdGradChi )[ numThirdOrderTensor * i + j ];

                }

            }

            // Set the Jacobians of the state variables
            std::copy( dStateVariableResidualsdD->begin( ), dStateVariableResidualsdD->end( ), std::begin( *dRdD.value ) + numConfigurationUnknowns * numConfigurationUnknowns );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );
            dRdT.zero( *getNumEquations( ) );

        }

        void residual::projectSuggestedX( std::vector< floatType > &trialX,
                                          const std::vector< floatType > &Xp ){
            /*!
             * Project the suggested unknown vector to the allowable space
             *
             * Called whenever hydra calls updateUnknownVector. It is assumed that the
             * initial value as suggested by `residual::suggestInitialIterationValues` is
             * in the allowable space.
             *
             * \param &trialX: The trial value of X
             * \param &Xp: The previously accepted value of X
             */

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            const unsigned int plasticConfigurationIndex = *getPlasticConfigurationIndex( );

            const unsigned int numConfigurations         = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns  = *hydra->getConfigurationUnknownCount( );

            const std::vector< unsigned int > *stateVariableIndices = getStateVariableIndices( );

            const floatVector dx = trialX - Xp;

            const unsigned int plasticDeformationStart = numConfigurationUnknowns * plasticConfigurationIndex;

            const unsigned int plasticDeformationStop  = numConfigurationUnknowns * ( plasticConfigurationIndex + 1 );

            floatVector deltaPlasticDeformations( dx.begin( ) + plasticDeformationStart,
                                                  dx.begin( ) + plasticDeformationStop );

            // Check the macro plastic deformation
            floatType norm = tardigradeVectorTools::l2norm( secondOrderTensor( deltaPlasticDeformations.begin( ),
                                                                               deltaPlasticDeformations.begin( ) + numSecondOrderTensor ) );

            if ( norm > *getMaxMacroPlasticDeltaNorm( ) ){

                for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                    trialX[ plasticDeformationStart + i ] = Xp[ i + plasticDeformationStart ] + dx[ i + plasticDeformationStart ] / norm;

                }

            }

            // Check the micro plastic deformation
            norm = tardigradeVectorTools::l2norm( secondOrderTensor( deltaPlasticDeformations.begin( ) + numSecondOrderTensor,
                                                                     deltaPlasticDeformations.begin( ) + 2 * numSecondOrderTensor ) );

            if ( norm > *getMaxMicroPlasticDeltaNorm( ) ){

                for ( unsigned int i = 0; i < numSecondOrderTensor; i++ ){

                    trialX[ plasticDeformationStart + i + numSecondOrderTensor ] = Xp[ plasticDeformationStart + i + numSecondOrderTensor ] + dx[ plasticDeformationStart + i + numSecondOrderTensor ] / norm;

                }

            }

            // Check the micro gradient plastic deformation
            norm = tardigradeVectorTools::l2norm( thirdOrderTensor( deltaPlasticDeformations.begin( ) + 2 * numSecondOrderTensor,
                                                                    deltaPlasticDeformations.begin( ) + 2 * numSecondOrderTensor + numThirdOrderTensor ) );

            if ( norm > *getMaxMicroGradientPlasticDeltaNorm( ) ){

                for ( unsigned int i = 0; i < numThirdOrderTensor; i++ ){

                    trialX[ plasticDeformationStart + i + 2 * numSecondOrderTensor ] = Xp[ plasticDeformationStart + i + 2 * numSecondOrderTensor ] + dx[ plasticDeformationStart + i + 2 * numSecondOrderTensor ] / norm;

                }

            }

            // Check the plastic multipliers and state variables to ensure positive values
            for ( auto v = stateVariableIndices->begin( ); v != stateVariableIndices->end( ); v++ ){

                if ( trialX[ numConfigurations * numConfigurationUnknowns + ( *v ) ] < 0 ){

                    trialX[ numConfigurations * numConfigurationUnknowns + ( *v ) ] = 0;

                }

            }

        }

    }

}
