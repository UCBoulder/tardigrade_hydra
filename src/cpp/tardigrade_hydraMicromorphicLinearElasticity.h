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
                    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( extractMaterialParameters( parameters, _Amatrix, _Bmatrix, _Cmatrix, _Dmatrix ) );

                }

                //! Return a reference to the A stiffness matrix
                const parameterVector *getAMatrix( ){ return &_Amatrix; }

                //! Return a reference to the B stiffness matrix
                const parameterVector *getBMatrix( ){ return &_Bmatrix; }

                //! Return a reference to the C stiffness matrix
                const parameterVector *getCMatrix( ){ return &_Cmatrix; }

                //! Return a reference to the D stiffness matrix
                const parameterVector *getDMatrix( ){ return &_Dmatrix; }

            protected:

                //! Set the current values of the deformation
                virtual void setDeformation( ){ setDeformation( false ); }

                //! Set the previous values of the deformation
                virtual void setPreviousDeformation( ){ setDeformation( true ); }

                virtual void setDeformation( const bool isPrevious );

                virtual void setRightCauchyGreen( );

                virtual void setPsi( );

                virtual void setGamma( );

                virtual void setPreviousRightCauchyGreen( );

                virtual void setPreviousPsi( );

                virtual void setPreviousGamma( );

                virtual void setReferenceStresses( const bool isPrevious );

                virtual void setPK2Stress( );

                virtual void setReferenceSymmetricMicroStress( );

                virtual void setReferenceHigherOrderStress( );

                virtual void setPreviousPK2Stress( );

                virtual void setPreviousReferenceSymmetricMicroStress( );

                virtual void setPreviousReferenceHigherOrderStress( );

                virtual void setStresses( const bool isPrevious );

                virtual void setCauchyStress( );

                virtual void setSymmetricMicroStress( );

                virtual void setHigherOrderStress( );

                virtual void setPreviousCauchyStress( );

                virtual void setPreviousSymmetricMicroStress( );

                virtual void setPreviousHigherOrderStress( );

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

                virtual void setdPK2dF( );

                virtual void setdPK2dFn( );

                virtual void setdPK2dChi( );

                virtual void setdPK2dChin( );

                virtual void setdPK2dGradChi( );

                virtual void setdPK2dGradChin( );

                virtual void setdSIGMAdF( );

                virtual void setdSIGMAdFn( );

                virtual void setdSIGMAdChi( );

                virtual void setdSIGMAdChin( );

                virtual void setdSIGMAdGradChi( );

                virtual void setdSIGMAdGradChin( );

                virtual void setdMdF( );

                virtual void setdMdFn( );

                virtual void setdMdChi( );

                virtual void setdMdChin( );

                virtual void setdMdGradChi( );

                virtual void setdMdGradChin( );

                virtual void setPreviousdPK2dF( );

                virtual void setPreviousdPK2dFn( );

                virtual void setPreviousdPK2dChi( );

                virtual void setPreviousdPK2dChin( );

                virtual void setPreviousdPK2dGradChi( );

                virtual void setPreviousdPK2dGradChin( );

                virtual void setPreviousdSIGMAdF( );

                virtual void setPreviousdSIGMAdFn( );

                virtual void setPreviousdSIGMAdChi( );

                virtual void setPreviousdSIGMAdChin( );

                virtual void setPreviousdSIGMAdGradChi( );

                virtual void setPreviousdSIGMAdGradChin( );

                virtual void setPreviousdMdF( );

                virtual void setPreviousdMdFn( );

                virtual void setPreviousdMdChi( );

                virtual void setPreviousdMdChin( );

                virtual void setPreviousdMdGradChi( );

                virtual void setPreviousdMdGradChin( );

                virtual void setReferenceStressJacobians( const bool isPrevious );

                virtual void setdCauchyStressdF( );

                virtual void setdCauchyStressdFn( );

                virtual void setdCauchyStressdChi( );

                virtual void setdCauchyStressdChin( );

                virtual void setdCauchyStressdGradChi( );

                virtual void setdCauchyStressdGradChin( );

                virtual void setdSymmetricMicroStressdF( );

                virtual void setdSymmetricMicroStressdFn( );

                virtual void setdSymmetricMicroStressdChi( );

                virtual void setdSymmetricMicroStressdChin( );

                virtual void setdSymmetricMicroStressdGradChi( );

                virtual void setdSymmetricMicroStressdGradChin( );

                virtual void setdHigherOrderStressdF( );

                virtual void setdHigherOrderStressdFn( );

                virtual void setdHigherOrderStressdChi( );

                virtual void setdHigherOrderStressdChin( );

                virtual void setdHigherOrderStressdGradChi( );

                virtual void setdHigherOrderStressdGradChin( );

                virtual void setPreviousdCauchyStressdF( );

                virtual void setPreviousdCauchyStressdFn( );

                virtual void setPreviousdCauchyStressdChi( );

                virtual void setPreviousdCauchyStressdChin( );

                virtual void setPreviousdCauchyStressdGradChi( );

                virtual void setPreviousdCauchyStressdGradChin( );

                virtual void setPreviousdSymmetricMicroStressdF( );

                virtual void setPreviousdSymmetricMicroStressdFn( );

                virtual void setPreviousdSymmetricMicroStressdChi( );

                virtual void setPreviousdSymmetricMicroStressdChin( );

                virtual void setPreviousdSymmetricMicroStressdGradChi( );

                virtual void setPreviousdSymmetricMicroStressdGradChin( );

                virtual void setPreviousdHigherOrderStressdF( );

                virtual void setPreviousdHigherOrderStressdFn( );

                virtual void setPreviousdHigherOrderStressdChi( );

                virtual void setPreviousdHigherOrderStressdChin( );

                virtual void setPreviousdHigherOrderStressdGradChi( );

                virtual void setPreviousdHigherOrderStressdGradChin( );

                virtual void setStressesJacobians( const bool isPrevious );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdD( ) override;

            private:
                // Friend classes
                friend class tardigradeHydra::micromorphicLinearElasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

                using tardigradeHydra::residualBaseMicromorphic::setResidual;

                using tardigradeHydra::residualBaseMicromorphic::setJacobian;

                using tardigradeHydra::residualBaseMicromorphic::setdRdF;

                using tardigradeHydra::residualBaseMicromorphic::setdRdT;

                using tardigradeHydra::residualBaseMicromorphic::setAdditionalDerivatives;

                using tardigradeHydra::residualBaseMicromorphic::setStress;

                using tardigradeHydra::residualBaseMicromorphic::setPreviousStress;

                parameterVector _Amatrix; //!< The A stiffness matrix

                parameterVector _Bmatrix; //!< The B stiffness matrix

                parameterVector _Cmatrix; //!< The C stiffness matrix

                parameterVector _Dmatrix; //!< The D stiffness matrix

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, rightCauchyGreen,                       floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dRightCauchyGreendF,                    floatMatrix, setdRightCauchyGreendF                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dRightCauchyGreendFn,                   floatMatrix, setdRightCauchyGreendFn                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, psi,                                    floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidF,                                 floatMatrix, setdPsidF                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidFn,                                floatMatrix, setdPsidFn                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidChi,                               floatMatrix, setdPsidChi                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidChin,                              floatMatrix, setdPsidChin                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, gamma,                                  floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadF,                               floatMatrix, setdGammadF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadChi,                             floatMatrix, setdGammadChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadGradChi,                         floatMatrix, setdGammadGradChi                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadFn,                              floatMatrix, setdGammadFn                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadChin,                            floatMatrix, setdGammadChin                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadGradChin,                        floatMatrix, setdGammadGradChin                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousRightCauchyGreen,               floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdRightCauchyGreendF,            floatMatrix, setPreviousdRightCauchyGreendF            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdRightCauchyGreendFn,           floatMatrix, setPreviousdRightCauchyGreendFn           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPsi,                            floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidF,                         floatMatrix, setPreviousdPsidF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidFn,                        floatMatrix, setPreviousdPsidFn                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidChi,                       floatMatrix, setPreviousdPsidChi                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidChin,                      floatMatrix, setPreviousdPsidChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousGamma,                          floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadF,                       floatMatrix, setPreviousdGammadF                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadChi,                     floatMatrix, setPreviousdGammadChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadGradChi,                 floatMatrix, setPreviousdGammadGradChi                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadFn,                      floatMatrix, setPreviousdGammadFn                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadChin,                    floatMatrix, setPreviousdGammadChin                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadGradChin,                floatMatrix, setPreviousdGammadGradChin                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, PK2Stress,                              floatVector, setPK2Stress                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dF,                                 floatMatrix, setdPK2dF                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dFn,                                floatMatrix, setdPK2dFn                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dChi,                               floatMatrix, setdPK2dChi                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dChin,                              floatMatrix, setdPK2dChin                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dGradChi,                           floatMatrix, setdPK2dGradChi                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dGradChin,                          floatMatrix, setdPK2dGradChin                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, referenceSymmetricMicroStress,          floatVector, setReferenceSymmetricMicroStress          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdF,                               floatMatrix, setdSIGMAdF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdFn,                              floatMatrix, setdSIGMAdFn                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdChi,                             floatMatrix, setdSIGMAdChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdChin,                            floatMatrix, setdSIGMAdChin                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdGradChi,                         floatMatrix, setdSIGMAdGradChi                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdGradChin,                        floatMatrix, setdSIGMAdGradChin                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, referenceHigherOrderStress,             floatVector, setReferenceHigherOrderStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdF,                                   floatMatrix, setdMdF                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdFn,                                  floatMatrix, setdMdFn                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdChi,                                 floatMatrix, setdMdChi                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdChin,                                floatMatrix, setdMdChin                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdGradChi,                             floatMatrix, setdMdGradChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdGradChin,                            floatMatrix, setdMdGradChin                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPK2Stress,                      floatVector, setPreviousPK2Stress                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dF,                         floatMatrix, setPreviousdPK2dF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dFn,                        floatMatrix, setPreviousdPK2dFn                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dChi,                       floatMatrix, setPreviousdPK2dChi                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dChin,                      floatMatrix, setPreviousdPK2dChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dGradChi,                   floatMatrix, setPreviousdPK2dGradChi                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dGradChin,                  floatMatrix, setPreviousdPK2dGradChin                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousReferenceSymmetricMicroStress,  floatVector, setPreviousReferenceSymmetricMicroStress  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdF,                       floatMatrix, setPreviousdSIGMAdF                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdFn,                      floatMatrix, setPreviousdSIGMAdFn                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdChi,                     floatMatrix, setPreviousdSIGMAdChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdChin,                    floatMatrix, setPreviousdSIGMAdChin                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdGradChi,                 floatMatrix, setPreviousdSIGMAdGradChi                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdGradChin,                floatMatrix, setPreviousdSIGMAdGradChin                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousReferenceHigherOrderStress,     floatVector, setPreviousReferenceHigherOrderStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdF,                           floatMatrix, setPreviousdMdF                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdFn,                          floatMatrix, setPreviousdMdFn                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdChi,                         floatMatrix, setPreviousdMdChi                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdChin,                        floatMatrix, setPreviousdMdChin                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdGradChi,                     floatMatrix, setPreviousdMdGradChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdGradChin,                    floatMatrix, setPreviousdMdGradChin                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, cauchyStress,                           floatVector, setCauchyStress                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdF,                        floatMatrix, setdCauchyStressdF                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdFn,                       floatMatrix, setdCauchyStressdFn                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdChi,                      floatMatrix, setdCauchyStressdChi                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdChin,                     floatMatrix, setdCauchyStressdChin                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdGradChi,                  floatMatrix, setdCauchyStressdGradChi                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdGradChin,                 floatMatrix, setdCauchyStressdGradChin                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, symmetricMicroStress,                   floatVector, setSymmetricMicroStress                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdF,                floatMatrix, setdSymmetricMicroStressdF                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdFn,               floatMatrix, setdSymmetricMicroStressdFn               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdChi,              floatMatrix, setdSymmetricMicroStressdChi              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdChin,             floatMatrix, setdSymmetricMicroStressdChin             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdGradChi,          floatMatrix, setdSymmetricMicroStressdGradChi          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdGradChin,         floatMatrix, setdSymmetricMicroStressdGradChin         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, higherOrderStress,                      floatVector, setHigherOrderStress                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdF,                   floatMatrix, setdHigherOrderStressdF                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdFn,                  floatMatrix, setdHigherOrderStressdFn                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdChi,                 floatMatrix, setdHigherOrderStressdChi                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdChin,                floatMatrix, setdHigherOrderStressdChin                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdGradChi,             floatMatrix, setdHigherOrderStressdGradChi             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdGradChin,            floatMatrix, setdHigherOrderStressdGradChin            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousCauchyStress,                   floatVector, setPreviousCauchyStress                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdF,                floatMatrix, setPreviousdCauchyStressdF                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdFn,               floatMatrix, setPreviousdCauchyStressdFn               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdChi,              floatMatrix, setPreviousdCauchyStressdChi              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdChin,             floatMatrix, setPreviousdCauchyStressdChin             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdGradChi,          floatMatrix, setPreviousdCauchyStressdGradChi          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdGradChin,         floatMatrix, setPreviousdCauchyStressdGradChin         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSymmetricMicroStress,           floatVector, setPreviousSymmetricMicroStress           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdF,        floatMatrix, setPreviousdSymmetricMicroStressdF        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdFn,       floatMatrix, setPreviousdSymmetricMicroStressdFn       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdChi,      floatMatrix, setPreviousdSymmetricMicroStressdChi      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdChin,     floatMatrix, setPreviousdSymmetricMicroStressdChin     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdGradChi,  floatMatrix, setPreviousdSymmetricMicroStressdGradChi  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdGradChin, floatMatrix, setPreviousdSymmetricMicroStressdGradChin )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHigherOrderStress,              floatVector, setPreviousHigherOrderStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdF,           floatMatrix, setPreviousdHigherOrderStressdF           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdFn,          floatMatrix, setPreviousdHigherOrderStressdFn          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdChi,         floatMatrix, setPreviousdHigherOrderStressdChi         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdChin,        floatMatrix, setPreviousdHigherOrderStressdChin        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdGradChi,     floatMatrix, setPreviousdHigherOrderStressdGradChi     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdGradChin,    floatMatrix, setPreviousdHigherOrderStressdGradChin    )

        };

    }

}

#endif
