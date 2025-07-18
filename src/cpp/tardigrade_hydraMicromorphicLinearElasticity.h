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

        void linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress );
    
        void linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress,
                               variableVector &dCauchyStressdF, variableVector &dCauchyStressdChi, variableVector &dCauchyStressdGradChi,
                               variableVector &dMicroStressdF, variableVector &dMicroStressdChi, variableVector &dMicroStressdGradChi,
                               variableVector &dHigherOrderStressdF, variableVector &dHigherOrderStressdChi,
                               variableVector &dHigherOrderStressdGradChi );

        void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress );
    
        void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableVector &dPK2StressdF, variableVector &dPK2StressdChi, variableVector &dPK2StressdGradChi,
                                        variableVector &dReferenceMicroStressdF, variableVector &dReferenceMicroStressdChi,
                                        variableVector &dReferenceMicroStressdGradChi, variableVector &dMdF, variableVector &dMdGradChi );

        void linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                           const variableVector &Gamma,
                                                           const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                           const parameterVector &D,
                                                           variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                           variableVector &referenceHigherOrderStress );
    
        void linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                           const variableVector &Gamma,
                                                           const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                           const parameterVector &D,
                                                           variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                           variableVector &referenceHigherOrderStress,
                                                           variableVector &dPK2StressdRCG, variableVector &dPK2StressdPsi,
                                                           variableVector &dPK2StressdGamma,
                                                           variableVector &dReferenceMicroStressdRCG,
                                                           variableVector &dReferenceMicroStressdPsi,
                                                           variableVector &dReferenceMicroStressdGamma,
                                                           variableVector &dMdGamma );

        void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress );
    
        void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress,
                                             variableVector &dCauchyStressdF, variableVector &dCauchyStressdPK2Stress,
                                             variableVector &dMicroStressdF, variableVector &dMicroStressdReferenceMicroStress,
                                             variableVector &dHigherOrderStressdF, variableVector &dHigherOrderStressdChi,
                                             variableVector &dHigherOrderStressdReferenceHigherOrderStress );

        void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma );
    
        void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                             variableVector &dCdF, variableVector &dPsidF, variableVector &dPsidChi,
                                             variableVector &dGammadF, variableVector &dGammadGradChi );

        void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1 );
    
        void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1,
                                            variableVector &dTerm1dGreenLagrangeStrain, variableVector &dTerm1dMicroStrain );
    
        void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &incCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2 );
    
        void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2, variableVector &dTerm2dGreenLagrangeStrain,
                                            variableVector &dTerm2dMicroStrain, variableVector &dTerm2dInvCPsi );

        void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress );
    
        void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress,
                                                    variableVector &dHigherOrderStressdGamma );

        void computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3 );
    
        void computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3,
                                            variableVector &dTerm3dInvCGamma, variableVector &dTerm3dReferenceHigherOrderStress );

        void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi );
    
        void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                                   variableVector &dInvRCGPsidRGG, variableVector &dInvRCGPsidPsi );
    
        void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma );
    
        void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                     variableVector &dInvRCGGammadRCG, variableVector &dInvRCGGammadGamma );

        void formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A );
    
        void formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                                 const parameterType &nu,  const parameterType &sigma, parameterVector &B );
    
        void formIsotropicC( const parameterVector &taus, parameterVector &C );
    
        void formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D );

        void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation );
    
        void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation,
                                                         variableVector &dFdGradU, variableVector &dChidPhi,
                                                         variableVector &dGradChidGradPhi );

        void extractMaterialParameters( const std::vector< double > &fparams,
                                            parameterVector &Amatrix, parameterVector &Bmatrix,
                                            parameterVector &Cmatrix, parameterVector &Dmatrix );

        /*!
         * The residual for a micromorphic linear elasticity constitutive equation
         */
        class residual : public tardigradeHydra::residualBaseMicromorphic {

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

            protected:

                using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

                using tardigradeHydra::residualBaseMicromorphic::setResidual;

                using tardigradeHydra::residualBaseMicromorphic::setJacobian;

                using tardigradeHydra::residualBaseMicromorphic::setdRdD;

                using tardigradeHydra::residualBaseMicromorphic::setdRdT;

                using tardigradeHydra::residualBaseMicromorphic::setAdditionalDerivatives;

                using tardigradeHydra::residualBaseMicromorphic::setStress;

                using tardigradeHydra::residualBaseMicromorphic::setPreviousStress;

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

                virtual void setdRdT( ) override;

                virtual void setStress( ) override;

                virtual void setPreviousStress( ) override;

            private:
                // Friend classes
                friend class tardigradeHydra::micromorphicLinearElasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                parameterVector _Amatrix; //!< The A stiffness matrix

                parameterVector _Bmatrix; //!< The B stiffness matrix

                parameterVector _Cmatrix; //!< The C stiffness matrix

                parameterVector _Dmatrix; //!< The D stiffness matrix

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, rightCauchyGreen,                       floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dRightCauchyGreendF,                    floatVector, setdRightCauchyGreendF                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dRightCauchyGreendFn,                   floatVector, setdRightCauchyGreendFn                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, psi,                                    floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidF,                                 floatVector, setdPsidF                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidFn,                                floatVector, setdPsidFn                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidChi,                               floatVector, setdPsidChi                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPsidChin,                              floatVector, setdPsidChin                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, gamma,                                  floatVector, setDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadF,                               floatVector, setdGammadF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadChi,                             floatVector, setdGammadChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadGradChi,                         floatVector, setdGammadGradChi                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadFn,                              floatVector, setdGammadFn                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadChin,                            floatVector, setdGammadChin                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGammadGradChin,                        floatVector, setdGammadGradChin                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousRightCauchyGreen,               floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdRightCauchyGreendF,            floatVector, setPreviousdRightCauchyGreendF            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdRightCauchyGreendFn,           floatVector, setPreviousdRightCauchyGreendFn           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPsi,                            floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidF,                         floatVector, setPreviousdPsidF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidFn,                        floatVector, setPreviousdPsidFn                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidChi,                       floatVector, setPreviousdPsidChi                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPsidChin,                      floatVector, setPreviousdPsidChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousGamma,                          floatVector, setPreviousDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadF,                       floatVector, setPreviousdGammadF                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadChi,                     floatVector, setPreviousdGammadChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadGradChi,                 floatVector, setPreviousdGammadGradChi                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadFn,                      floatVector, setPreviousdGammadFn                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadChin,                    floatVector, setPreviousdGammadChin                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdGammadGradChin,                floatVector, setPreviousdGammadGradChin                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, PK2Stress,                              floatVector, setPK2Stress                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dF,                                 floatVector, setdPK2dF                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dFn,                                floatVector, setdPK2dFn                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dChi,                               floatVector, setdPK2dChi                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dChin,                              floatVector, setdPK2dChin                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dGradChi,                           floatVector, setdPK2dGradChi                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2dGradChin,                          floatVector, setdPK2dGradChin                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, referenceSymmetricMicroStress,          floatVector, setReferenceSymmetricMicroStress          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdF,                               floatVector, setdSIGMAdF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdFn,                              floatVector, setdSIGMAdFn                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdChi,                             floatVector, setdSIGMAdChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdChin,                            floatVector, setdSIGMAdChin                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdGradChi,                         floatVector, setdSIGMAdGradChi                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSIGMAdGradChin,                        floatVector, setdSIGMAdGradChin                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, referenceHigherOrderStress,             floatVector, setReferenceHigherOrderStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdF,                                   floatVector, setdMdF                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdFn,                                  floatVector, setdMdFn                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdChi,                                 floatVector, setdMdChi                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdChin,                                floatVector, setdMdChin                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdGradChi,                             floatVector, setdMdGradChi                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMdGradChin,                            floatVector, setdMdGradChin                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPK2Stress,                      floatVector, setPreviousPK2Stress                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dF,                         floatVector, setPreviousdPK2dF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dFn,                        floatVector, setPreviousdPK2dFn                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dChi,                       floatVector, setPreviousdPK2dChi                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dChin,                      floatVector, setPreviousdPK2dChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dGradChi,                   floatVector, setPreviousdPK2dGradChi                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2dGradChin,                  floatVector, setPreviousdPK2dGradChin                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousReferenceSymmetricMicroStress,  floatVector, setPreviousReferenceSymmetricMicroStress  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdF,                       floatVector, setPreviousdSIGMAdF                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdFn,                      floatVector, setPreviousdSIGMAdFn                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdChi,                     floatVector, setPreviousdSIGMAdChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdChin,                    floatVector, setPreviousdSIGMAdChin                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdGradChi,                 floatVector, setPreviousdSIGMAdGradChi                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSIGMAdGradChin,                floatVector, setPreviousdSIGMAdGradChin                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousReferenceHigherOrderStress,     floatVector, setPreviousReferenceHigherOrderStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdF,                           floatVector, setPreviousdMdF                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdFn,                          floatVector, setPreviousdMdFn                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdChi,                         floatVector, setPreviousdMdChi                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdChin,                        floatVector, setPreviousdMdChin                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdGradChi,                     floatVector, setPreviousdMdGradChi                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMdGradChin,                    floatVector, setPreviousdMdGradChin                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, cauchyStress,                           floatVector, setCauchyStress                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdF,                        floatVector, setdCauchyStressdF                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdFn,                       floatVector, setdCauchyStressdFn                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdChi,                      floatVector, setdCauchyStressdChi                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdChin,                     floatVector, setdCauchyStressdChin                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdGradChi,                  floatVector, setdCauchyStressdGradChi                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdGradChin,                 floatVector, setdCauchyStressdGradChin                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, symmetricMicroStress,                   floatVector, setSymmetricMicroStress                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdF,                floatVector, setdSymmetricMicroStressdF                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdFn,               floatVector, setdSymmetricMicroStressdFn               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdChi,              floatVector, setdSymmetricMicroStressdChi              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdChin,             floatVector, setdSymmetricMicroStressdChin             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdGradChi,          floatVector, setdSymmetricMicroStressdGradChi          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroStressdGradChin,         floatVector, setdSymmetricMicroStressdGradChin         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, higherOrderStress,                      floatVector, setHigherOrderStress                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdF,                   floatVector, setdHigherOrderStressdF                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdFn,                  floatVector, setdHigherOrderStressdFn                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdChi,                 floatVector, setdHigherOrderStressdChi                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdChin,                floatVector, setdHigherOrderStressdChin                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdGradChi,             floatVector, setdHigherOrderStressdGradChi             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderStressdGradChin,            floatVector, setdHigherOrderStressdGradChin            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousCauchyStress,                   floatVector, setPreviousCauchyStress                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdF,                floatVector, setPreviousdCauchyStressdF                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdFn,               floatVector, setPreviousdCauchyStressdFn               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdChi,              floatVector, setPreviousdCauchyStressdChi              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdChin,             floatVector, setPreviousdCauchyStressdChin             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdGradChi,          floatVector, setPreviousdCauchyStressdGradChi          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdGradChin,         floatVector, setPreviousdCauchyStressdGradChin         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSymmetricMicroStress,           floatVector, setPreviousSymmetricMicroStress           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdF,        floatVector, setPreviousdSymmetricMicroStressdF        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdFn,       floatVector, setPreviousdSymmetricMicroStressdFn       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdChi,      floatVector, setPreviousdSymmetricMicroStressdChi      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdChin,     floatVector, setPreviousdSymmetricMicroStressdChin     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdGradChi,  floatVector, setPreviousdSymmetricMicroStressdGradChi  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroStressdGradChin, floatVector, setPreviousdSymmetricMicroStressdGradChin )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHigherOrderStress,              floatVector, setPreviousHigherOrderStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdF,           floatVector, setPreviousdHigherOrderStressdF           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdFn,          floatVector, setPreviousdHigherOrderStressdFn          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdChi,         floatVector, setPreviousdHigherOrderStressdChi         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdChin,        floatVector, setPreviousdHigherOrderStressdChin        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdGradChi,     floatVector, setPreviousdHigherOrderStressdGradChi     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderStressdGradChin,    floatVector, setPreviousdHigherOrderStressdGradChin    )

        };

    }

}

#endif
