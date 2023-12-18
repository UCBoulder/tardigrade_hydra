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
             * A_{KLMN} = \lambda \delta_{KL} \delta_{MN} + \mu \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )
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
             * B_{KLMN} = ( eta - tau ) \delta_{KL} \delta_{MN} + \kappa \delta_{KM} \delta_{LN} + \nu \delta_{KN} \delta_{LM}
             *          - \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right)
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
             * C_{KLMNPQ} = \tau_1 \left( \delta_{KL} \delta_{MN} \delta_{PQ} + \delta_{KQ} \delta_{LM} \delta_{NP} \right) 
             *            + \tau_2 \left( \delta_{KL} \delta_{MP} \delta_{NQ} + \delta_{KM} \delta_{LQ} \delta_{NP} \right)
             *            + \tau_3 \delta_{KL} \delta_{MQ} \delta_{NP}
             *            + \tau_4 \delta_{KN} \delta_{LM} \delta_{PQ}
             *            + \tau_5 \left( \delta_{KM} \delta_{LN} \delta_{PQ} + \delta_{KP} \delta_{LM} \delta_{NQ} )
             *            + \tau_6 \delta_{KM} \delta_{LP} \delta_{NQ}
             *            + \tau_7 \delta_{KN} \delta_{LP} \delta_{MQ}
             *            + \tau_8 \left( \delta_{KP} \delta_{LQ} \delta_{MN} + \delta_{KQ} \delta_{LN} \delta_{MP} )
             *            + \tau_9 \delta_{KN} \delta_{LQ} \delta_{MP}
             *            + \tau_{10} \delta_{KP} \delta_{LN} \delta_{MQ}
             *            + \tau_{11} \delta_{KQ} \delta_{LP} \delta_{MN}
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
             * D_{KLMN} = \tau \delta_{KL} \delta_{MN} + \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )
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

    }

}
