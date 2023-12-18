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
             * B_{KLMN} &= ( eta - tau ) \delta_{KL} \delta_{MN} + \kappa \delta_{KM} \delta_{LN} + \nu \delta_{KN} \delta_{LM}
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
             * C_{KLMNPQ} &= \tau_1 \left( \delta_{KL} \delta_{MN} \delta_{PQ} + \delta_{KQ} \delta_{LM} \delta_{NP} \right) 
             *            &+ \tau_2 \left( \delta_{KL} \delta_{MP} \delta_{NQ} + \delta_{KM} \delta_{LQ} \delta_{NP} \right)
             *            &+ \tau_3 \delta_{KL} \delta_{MQ} \delta_{NP}
             *            &+ \tau_4 \delta_{KN} \delta_{LM} \delta_{PQ}
             *            &+ \tau_5 \left( \delta_{KM} \delta_{LN} \delta_{PQ} + \delta_{KP} \delta_{LM} \delta_{NQ} )
             *            &+ \tau_6 \delta_{KM} \delta_{LP} \delta_{NQ}
             *            &+ \tau_7 \delta_{KN} \delta_{LP} \delta_{MQ}
             *            &+ \tau_8 \left( \delta_{KP} \delta_{LQ} \delta_{MN} + \delta_{KQ} \delta_{LN} \delta_{MP} )
             *            &+ \tau_9 \delta_{KN} \delta_{LQ} \delta_{MP}
             *            &+ \tau_{10} \delta_{KP} \delta_{LN} \delta_{MQ}
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

    }

}
