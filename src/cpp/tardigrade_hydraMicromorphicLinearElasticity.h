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
#include<tardigrade_hydraMicromorphic.h>

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

        errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                                   const variableVector &gradientMicroDeformation,
                                   const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                   const parameterVector &D,
                                   variableVector &cauchyStress, variableVector &microStress,
                                   variableVector &higherOrderStress );
    
        errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                                   const variableVector &gradientMicroDeformation,
                                   const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                   const parameterVector &D,
                                   variableVector &cauchyStress, variableVector &microStress,
                                   variableVector &higherOrderStress,
                                   variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdChi, variableMatrix &dCauchyStressdGradChi,
                                   variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdChi, variableMatrix &dMicroStressdGradChi,
                                   variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                   variableMatrix &dHigherOrderStressdGradChi );

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

        class residual : public tardigradeHydra::residualBaseMicromorphic {
            /*!
             * The residual for a micromorphic linear elasticity constitutive equation
             */

            public:

                residual( hydraBaseMicromorphic *_hydra, const unsigned int &_numEquations, const floatVector &parameters ) : tardigradeHydra::residualBaseMicromorphic( _hydra, _numEquations ){
                    /*!
                     * The main initialization constructor for the linear elastic residual
                     * 
                     * \param *_hydra: A pointer to the containing hydra class
                     * \param &_numEquations: The number of equations the residual defines
                     * \param &parameters: The parameter vector
                     */

                    // Form the stiffness matrices
                    TARDIGRADE_ERROR_TOOLS_CATCH( extractMaterialParameters( parameters, _Amatrix, _Bmatrix, _Cmatrix, _Dmatrix ) );

                }

                //! Return a reference to the A stiffness matrix
                const parameterVector *getAMatrix( ){ return &_Amatrix; }

                //! Return a reference to the B stiffness matrix
                const parameterVector *getBMatrix( ){ return &_Bmatrix; }

                //! Return a reference to the C stiffness matrix
                const parameterVector *getCMatrix( ){ return &_Cmatrix; }

                //! Return a reference to the D stiffness matrix
                const parameterVector *getDMatrix( ){ return &_Dmatrix; }

                const variableVector *getRightCauchyGreen( );

                const variableVector *getPsi( );

                const variableVector *getGamma( );

                const variableVector *getPreviousRightCauchyGreen( );

                const variableVector *getPreviousPsi( );

                const variableVector *getPreviousGamma( );

                const variableMatrix *getdRightCauchyGreendF( );

                const variableMatrix *getdRightCauchyGreendFn( );

                const variableMatrix *getdPsidF( );

                const variableMatrix *getdPsidFn( );

                const variableMatrix *getdPsidChi( );

                const variableMatrix *getdPsidChin( );

                const variableMatrix *getdGammadF( );

                const variableMatrix *getdGammadFn( );

                const variableMatrix *getdGammadChi( );

                const variableMatrix *getdGammadChin( );

                const variableMatrix *getdGammadGradChi( );

                const variableMatrix *getdGammadGradChin( );

                const variableMatrix *getPreviousdRightCauchyGreendF( );

                const variableMatrix *getPreviousdRightCauchyGreendFn( );

                const variableMatrix *getPreviousdPsidF( );

                const variableMatrix *getPreviousdPsidFn( );

                const variableMatrix *getPreviousdPsidChi( );

                const variableMatrix *getPreviousdPsidChin( );

                const variableMatrix *getPreviousdGammadF( );

                const variableMatrix *getPreviousdGammadFn( );

                const variableMatrix *getPreviousdGammadChi( );

                const variableMatrix *getPreviousdGammadChin( );

                const variableMatrix *getPreviousdGammadGradChi( );

                const variableMatrix *getPreviousdGammadGradChin( );

            protected:

                void setRightCauchyGreen( const variableVector &rightCauchyGreen );

                void setPsi( const variableVector &psi );

                void setGamma( const variableVector &gamma );

                void setPreviousRightCauchyGreen( const variableVector &previousRightCauchyGreen );

                void setPreviousPsi( const variableVector &previousPsi );

                void setPreviousGamma( const variableVector &previousGamma );

                void setdRightCauchyGreendF( const floatMatrix &value );

                void setdRightCauchyGreendFn( const floatMatrix &value );

                void setdPsidF( const floatMatrix &value );

                void setdPsidFn( const floatMatrix &value );

                void setdPsidChi( const floatMatrix &value );

                void setdPsidChin( const floatMatrix &value );

                void setdGammadF( const floatMatrix &value );

                void setdGammadFn( const floatMatrix &value );

                void setdGammadChi( const floatMatrix &value );

                void setdGammadChin( const floatMatrix &value );

                void setdGammadGradChi( const floatMatrix &value );

                void setdGammadGradChin( const floatMatrix &value );

                void setPreviousdRightCauchyGreendF( const floatMatrix &value );

                void setPreviousdRightCauchyGreendFn( const floatMatrix &value );

                void setPreviousdPsidF( const floatMatrix &value );

                void setPreviousdPsidFn( const floatMatrix &value );

                void setPreviousdPsidChi( const floatMatrix &value );

                void setPreviousdPsidChin( const floatMatrix &value );

                void setPreviousdGammadF( const floatMatrix &value );

                void setPreviousdGammadFn( const floatMatrix &value );

                void setPreviousdGammadChi( const floatMatrix &value );

                void setPreviousdGammadChin( const floatMatrix &value );

                void setPreviousdGammadGradChi( const floatMatrix &value );

                void setPreviousdGammadGradChin( const floatMatrix &value );

                virtual void setDeformation( const bool isPrevious );

                virtual void setRightCauchyGreen( );

                virtual void setPsi( );

                virtual void setGamma( );

                virtual void setPreviousRightCauchyGreen( );

                virtual void setPreviousPsi( );

                virtual void setPreviousGamma( );

                virtual void setDeformationJacobians( const bool isPrevious );

                virtual void setdRightCauchyGreendF( );

                virtual void setdRightCauchyGreendFn( );

                virtual void setdPsidF( );

                virtual void setdPsidFn( );

                virtual void setdPsidChi( );

                virtual void setdPsidChin( );

                virtual void setdGammadF( );

                virtual void setdGammadFn( );

                virtual void setdGammadChi( );

                virtual void setdGammadChin( );

                virtual void setdGammadGradChi( );

                virtual void setdGammadGradChin( );

                virtual void setPreviousdRightCauchyGreendF( );

                virtual void setPreviousdRightCauchyGreendFn( );

                virtual void setPreviousdPsidF( );

                virtual void setPreviousdPsidFn( );

                virtual void setPreviousdPsidChi( );

                virtual void setPreviousdPsidChin( );

                virtual void setPreviousdGammadF( );

                virtual void setPreviousdGammadFn( );

                virtual void setPreviousdGammadChi( );

                virtual void setPreviousdGammadChin( );

                virtual void setPreviousdGammadGradChi( );

                virtual void setPreviousdGammadGradChin( );

            private:

                parameterVector _Amatrix; //!< The A stiffness matrix

                parameterVector _Bmatrix; //!< The B stiffness matrix

                parameterVector _Cmatrix; //!< The C stiffness matrix

                parameterVector _Dmatrix; //!< The D stiffness matrix

                tardigradeHydra::dataStorage< variableVector > _rightCauchyGreen; //!< The current right Cauchy-Green deformation tensor

                tardigradeHydra::dataStorage< variableMatrix > _dRightCauchyGreendF; //!< The Jacobian of the right Cauchy-Green deformation tensor w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _dRightCauchyGreendFn; //!< The Jacobian of the right Cauchy-Green deformation tensor w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableVector > _psi; //!< The current micro-deformation tensor Psi

                tardigradeHydra::dataStorage< variableMatrix > _dPsidF; //!< The Jacobian of the micro deformation measure Psi w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _dPsidFn; //!< The Jacobian of the micro deformation measure Psi w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableMatrix > _dPsidChi; //!< The Jacobian of the micro deformation measure Psi w.r.t. the total micro deformation

                tardigradeHydra::dataStorage< variableMatrix > _dPsidChin; //!< The Jacobian of the micro deformation measure Psi w.r.t. the remaining sub micro deformations

                tardigradeHydra::dataStorage< variableVector > _gamma; //!< The current gradient micro-deformation tensor Gamma

                tardigradeHydra::dataStorage< variableMatrix > _dGammadF; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _dGammadChi; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the total micro-deformation

                tardigradeHydra::dataStorage< variableMatrix > _dGammadGradChi; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the reference spatial gradient of the total micro deformation

                tardigradeHydra::dataStorage< variableMatrix > _dGammadFn; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableMatrix > _dGammadChin; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the remaining sub micro-deformation

                tardigradeHydra::dataStorage< variableMatrix > _dGammadGradChin; //!< The Jacobian of the micro deformation measure Gamma w.r.t. the local reference spatial gradients of the remaining sub-micro deformations

                tardigradeHydra::dataStorage< variableVector > _previousRightCauchyGreen; //!< The previous right Cauchy-Green deformation tensor

                tardigradeHydra::dataStorage< variableMatrix > _previousdRightCauchyGreendF; //!< The Jacobian of the previous right Cauchy-Green deformation tensor w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _previousdRightCauchyGreendFn; //!< The Jacobian of the previous right Cauchy-Green deformation tensor w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableVector > _previousPsi; //!< The previous micro-deformation tensor Psi

                tardigradeHydra::dataStorage< variableMatrix > _previousdPsidF; //!< The Jacobian of the previous micro deformation measure Psi w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _previousdPsidFn; //!< The Jacobian of the previous micro deformation measure Psi w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableMatrix > _previousdPsidChi; //!< The Jacobian of the previous micro deformation measure Psi w.r.t. the total micro deformation

                tardigradeHydra::dataStorage< variableMatrix > _previousdPsidChin; //!< The Jacobian of the previous micro deformation measure Psi w.r.t. the remaining sub micro deformations

                tardigradeHydra::dataStorage< variableVector > _previousGamma; //!< The previous gradient micro-deformation tensor Gamma

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadF; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the total deformation gradient

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadChi; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the total micro-deformation

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadGradChi; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the reference spatial gradient of the total micro deformation

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadFn; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the remaining sub-deformation gradients

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadChin; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the remaining sub micro-deformation

                tardigradeHydra::dataStorage< variableMatrix > _previousdGammadGradChin; //!< The Jacobian of the previous micro deformation measure Gamma w.r.t. the local reference spatial gradients of the remaining sub-micro deformations
        };

    }

}

#endif
