/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicLinearElasticity.h
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework. Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_MICROMORPHIC_LINEAR_ELASTICITY_H
#define TARDIGRADE_HYDRA_MICROMORPHIC_LINEAR_ELASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace micromorphicLinearElasticity{

        // forward class definitions
        namespace unit_test{
            class residualTester;
        }

        constexpr const char* str_end(const char *str) {
            /*! Recursively search string for last character
             * \param *str: pointer to string START of UNIX path like string
             * \return *str: pointer to last character in string
             */
            return *str ? str_end(str + 1) : str;
        }
        constexpr bool str_slant(const char *str) {
            /*! Recursively search string for leftmost UNIX path separator from the left
             * \param *str: pointer to string START of UNIX path like string
             * \return bool: True if string contains UNIX path separator. Else false.
             */
            return *str == '/' ? true : (*str ? str_slant(str + 1) : false);
        }
        constexpr const char* r_slant(const char* str) {
            /*! Recursively search string for rightmost UNIX path separator from the right
             * \param *str: pointer to string END of UNIX path like string
             * \return *str: pointer to start of base name
             */
            return *str == '/' ? (str + 1) : r_slant(str - 1);
        }
        constexpr const char* file_name(const char* str) {
            /*! Return the current file name with extension at compile time
             * \param *str: pointer to string START of UNIX path like string
             * \return str: file base name
             */
            return str_slant(str) ? r_slant(str_end(str)) : str;
        }
        //Return filename for constructing debugging messages
        //https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
        const std::string __BASENAME__ = file_name(__FILE__);  //!< The base filename which will be parsed
        const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of(".")); //!< The parsed filename for error handling

        typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
        typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
        typedef double floatType; //!< Define the float values type.
        typedef std::vector< floatType > floatVector; //!< Define a vector of floats
        typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

        typedef floatType variableType; //!< Define the variable values type.
        typedef std::vector< variableType > variableVector; //!< Define a vector of variables
        typedef std::vector< std::vector< variableType > > variableMatrix; //!< Define a matrix of variables

        typedef double parameterType; //!< Define the parameter values type.
        typedef std::vector< parameterType > parameterVector; //!< Define a vector of parameters

        typedef double constantType; //!< Define the constant values type.
        typedef std::vector< constantType > constantVector; //!< Define a vector of constants
        typedef std::vector< std::vector< constantType > > constantMatrix; //!< Define a matrix of constants

        errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                            const variableVector &gradientMicroDeformation,
                                            const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                            const parameterVector &D,
                                            variableVector &PK2Stress, variableVector &referenceMicroStress,
                                            variableVector &referenceHigherOrderStress );
    
        errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                            const variableVector &gradientMicroDeformation,
                                            const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                            const parameterVector &D,
                                            variableVector &PK2Stress, variableVector &referenceMicroStress,
                                            variableVector &referenceHigherOrderStress,
                                            variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdChi, variableMatrix &dPK2StressdGradChi,
                                            variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdChi,
                                            variableMatrix &dReferenceMicroStressdGradChi, variableMatrix &dMdF, variableMatrix &dMdGradChi );

        errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                           const variableVector &Gamma,
                                                           const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                           const parameterVector &D,
                                                           variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                           variableVector &referenceHigherOrderStress );
    
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
                                                           variableMatrix &dMdGamma );

        errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress );
    
        errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress,
                                             variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                             variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                             variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                             variableMatrix &dHigherOrderStressdReferenceHigherOrderStress );

        errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma );
    
        errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                             variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidChi,
                                             variableMatrix &dGammadF, variableMatrix &dGammadGradChi );

        errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1 );
    
        errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1,
                                            variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain );
    
        errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &incCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2 );
    
        errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                            variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi );

        errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress );
    
        errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress,
                                                    variableMatrix &dHigherOrderStressdGamma );

        errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3 );
    
        errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3,
                                            variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress );

        errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi );
    
        errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                                   variableMatrix &dInvRCGPsidRGG, variableMatrix &dInvRCGPsidPsi );
    
        errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma );
    
        errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                     variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma );

        errorOut formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A );
    
        errorOut formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                                 const parameterType &nu,  const parameterType &sigma, parameterVector &B );
    
        errorOut formIsotropicC( const parameterVector &taus, parameterVector &C );
    
        errorOut formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D );

        errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation );
    
        errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation,
                                                         variableMatrix &dFdGradU, variableMatrix &dChidPhi,
                                                         variableMatrix &dGradChidGradPhi );

        errorOut extractMaterialParameters( const std::vector< double > &fparams,
                                            parameterVector &Amatrix, parameterVector &Bmatrix,
                                            parameterVector &Cmatrix, parameterVector &Dmatrix );

    }

}

#endif
