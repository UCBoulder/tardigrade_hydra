/**
  ******************************************************************************
  * \file tardigrade_hydraLinearViscoelasticity.h
  ******************************************************************************
  * An implementation of linear viscoelasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_LINEAR_VISCOELASTICITY_H
#define TARDIGRADE_HYDRA_LINEAR_VISCOELASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydraLinearElasticity.h>

namespace tardigradeHydra{

    namespace linearViscoelasticity{

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
        const std::string __BASENAME__ = file_name(__FILE__); //!< The base filename which will be parsed
    const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of(".")); //!< The parsed filename for error handling
    
        typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
        typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
        typedef double floatType; //!< Define the float values type.
        typedef std::vector< floatType > floatVector; //!< Define a vector of floats
        typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats 
    
        /*!
         * A residual class for a linear-viscoelastic material model where the stress is
         * computed in the reference configuration and pushed forward to the current
         * configuration.
         */
        class residual : public tardigradeHydra::linearElasticity::residual {
        
            public:

                /*!
                 * The main residual function for the linear viscoelastic residual
                 * 
                 * \param *hydra: The containing hydraBase object
                 * \param &numEquations: The number of equations the residual is responsible for
                 * \param &parameters: The parameter vector
                 * \param &viscoelasticISVLowerIndex: The lower index of the viscoelastic ISVs in the hydra::_additionalStateVariables
                 * \param &viscoelasticISVUpperIndex: The upper index (not included) of the viscoelastic ISVs in the hydra::_additionalStateVariables
                 * \param &integrationAlpha: The integration alpha parameter (0. is implicit, 1 is explicit)
                 */
                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const floatVector &parameters,
                          const unsigned int viscoelasticISVLowerIndex,
                          const unsigned int viscoelasticISVUpperIndex,
                          const floatType integrationAlpha=0. ) : tardigradeHydra::linearElasticity::residual( hydra, numEquations ), _viscoelasticISVLowerIndex( viscoelasticISVLowerIndex ), _viscoelasticISVUpperIndex( viscoelasticISVUpperIndex ), _integrationAlpha( integrationAlpha ){
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameterVector( parameters ) );
    
                }

                //! Get the lower index of the viscoelastic ISVs from the non-nonlinear solve state variable vector
                const unsigned int* getViscoelasticISVLowerIndex( ){ return &_viscoelasticISVLowerIndex; }

                //! Get the upper (but not including) index of the viscoelastic ISVs from the non-nonlinear solve state variable vector
                const unsigned int* getViscoelasticISVUpperIndex( ){ return &_viscoelasticISVUpperIndex; }

                //! Get the integration alpha parameter (0 is implicit, 1 is explicit)
                const floatType* getIntegrationAlpha( ){ return &_integrationAlpha; }

                //! Get the number of volumetric viscous terms
                const unsigned int* getNumVolumetricViscousTerms( ){ return &_numVolumetricViscousTerms; }

                //! Get the number of isochoric viscous terms
                const unsigned int* getNumIsochoricViscousTerms( ){ return &_numIsochoricViscousTerms; }

                //! Get the infinite bulk modulus
                const floatType* getKinf( ){ return &_Kinf; }

                //! Get the infinite shear modulus
                const floatType* getGinf( ){ return &_Ginf; }

                //! Get the bulk moduli
                const floatVector* getVolumetricModuli( ){ return &_Ks; }

                //! Get the isochoric moduli
                const floatVector* getIsochoricModuli( ){ return &_Gs; }

                //! Get the time constants for the volumetric moduli
                const floatVector* getVolumetricTaus( ){ return &_volumetricTaus; }

                //! Get the time constants for the isochoric moduli
                const floatVector* getIsochoricTaus( ){ return &_isochoricTaus; }

                virtual void decomposeDeformation( const floatVector &F, floatType &J, floatVector &Fhat );

                //! Get a pointer to the number of volumetric temperature parameters
                const floatVector* getVolumetricTemperatureParameters( ){ return &_volumetricTemperatureParameters; }

                //! Get a pointer to the number of isochoric temperature parameters
                const floatVector* getIsochoricTemperatureParameters( ){ return &_isochoricTemperatureParameters; }

                virtual void decomposeElasticDeformation( );

                virtual void decomposePreviousElasticDeformation( );

                virtual void setdJedFe( );

                virtual void setdFehatdFe( );

                void setNumStateVariables( const unsigned int numStateVariables );

                //! Get a pointer to the number of state variables
                unsigned int *getNumStateVariables( ){ return &_numStateVariables; }

                virtual void decomposeStateVariableVector( floatVector &volumetricStateVariables, floatVector &isochoricStateVariables );

                virtual floatType computeRateMultiplier( const floatVector &variables, const floatVector &parameters );

                virtual floatVector computedRateMultiplierdVariables( const floatVector &variables, const floatVector &parameters );

                void setVolumetricTemperatureParameters( const floatVector &parameters );

                void setIsochoricTemperatureParameters( const floatVector &parameters );

                virtual floatVector getVolumetricViscoelasticParameters( );

                virtual floatVector getIsochoricViscoelasticParameters( );

                void setdRdT( const floatVector &dRdT );

            protected:

                virtual void setNumVolumetricViscousTerms( const unsigned int &num );

                virtual void setNumIsochoricViscousTerms( const unsigned int &num );

                virtual void setKinf( const floatType &Kinf );

                virtual void setGinf( const floatType &Ginf );

                virtual void setVolumetricModuli( const floatVector &Ks );

                virtual void setIsochoricModuli( const floatVector &Gs );
 
                virtual void setVolumetricTaus( const floatVector &taus );

                virtual void setIsochoricTaus( const floatVector &taus );

                virtual void setPK2Stress( ) override;
    
                virtual void decomposeParameterVector( const floatVector &parameters );

                virtual void setVolumetricRateMultiplier( );

                virtual void setPreviousVolumetricRateMultiplier( );

                virtual void setIsochoricRateMultiplier( );

                virtual void setPreviousIsochoricRateMultiplier( );

                virtual void setdVolumetricRateMultiplierdT( );

                virtual void setdPreviousVolumetricRateMultiplierdPreviousT( );

                virtual void setdIsochoricRateMultiplierdT( );

                virtual void setdPreviousIsochoricRateMultiplierdPreviousT( );

                virtual void setPK2MeanStress( );

                virtual void setPK2IsochoricStress( );

                virtual void setdPK2MeanStressdFe( );

                virtual void setdPK2IsochoricStressdFe( );

                virtual void setdPK2MeanStressdT( );

                virtual void setdPK2IsochoricStressdT( );

                virtual void setdPK2StressdFe( ) override;

                virtual void setdPK2StressdT( );

                virtual void setdCauchyStressdT( );

                virtual void setdRdT( ) override;

                virtual void setUpdatedVolumetricViscoelasticStateVariables( );

                virtual void setUpdatedIsochoricViscoelasticStateVariables( );

                virtual void setCurrentAdditionalStateVariables( );

            private:

                using tardigradeHydra::linearElasticity::residual::residual;

                const unsigned int _viscoelasticISVLowerIndex; //!< The lower index of the viscoelastic ISVs

                const unsigned int _viscoelasticISVUpperIndex; //!< The not-included upper index of the viscoelastic ISVs

                const floatType _integrationAlpha; //!< The parameter for implicit (0) vs explicit (1) integration. Defaults to 0.

                unsigned int _numVolumetricViscousTerms; //!< The number of volumetric viscous terms

                unsigned int _numIsochoricViscousTerms; //!< The number of isochoric viscous terms

                unsigned int _numStateVariables; //!< The number of state variables required

                floatType _Kinf; //!< The infinite bulk modulus

                floatType _Ginf; //!< The infinite shear modulus

                floatVector _Ks; //!< The infinite bulk moduli

                floatVector _Gs; //!< The infinite shear moduli

                floatVector _volumetricTaus; //!< The volumetric time constants

                floatVector _isochoricTaus; //!< The isochoric time constants

                floatVector _volumetricTemperatureParameters; //!< The temperature parameters for the volumetric viscous elements

                floatVector _isochoricTemperatureParameters; //!< The temperature parameters for the isochoric viscous elements

                // Friend classes
                friend class tardigradeHydra::linearViscoelasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!
        
                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, Je,                                          floatType,   decomposeElasticDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, Fehat,                                       floatVector, decomposeElasticDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousJe,                                  floatType,   decomposePreviousElasticDeformation            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousFehat,                               floatVector, decomposePreviousElasticDeformation            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dJedFe,                                      floatVector, setdJedFe                                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFehatdFe,                                   floatMatrix, setdFehatdFe                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, volumetricRateMultiplier,                    floatType,   setVolumetricRateMultiplier                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousVolumetricRateMultiplier,            floatType,   setPreviousVolumetricRateMultiplier            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, isochoricRateMultiplier,                     floatType,   setIsochoricRateMultiplier                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousIsochoricRateMultiplier,             floatType,   setPreviousIsochoricRateMultiplier             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVolumetricRateMultiplierdT,                 floatType,   setdVolumetricRateMultiplierdT                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVolumetricRateMultiplierdPreviousT, floatType,   setdPreviousVolumetricRateMultiplierdPreviousT )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dIsochoricRateMultiplierdT,                  floatType,   setdIsochoricRateMultiplierdT                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousIsochoricRateMultiplierdPreviousT,  floatType,   setdPreviousIsochoricRateMultiplierdPreviousT  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, volumetricViscoelasticStateVariables,        floatVector, setUpdatedVolumetricViscoelasticStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, isochoricViscoelasticStateVariables,         floatVector, setUpdatedIsochoricViscoelasticStateVariables  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, PK2MeanStress,                               floatType,   setPK2MeanStress                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, PK2IsochoricStress,                          floatVector, setPK2IsochoricStress                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2MeanStressdT,                            floatType,   setdPK2MeanStressdT                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2MeanStressdFe,                           floatVector, setdPK2MeanStressdFe                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2IsochoricStressdT,                       floatVector, setdPK2IsochoricStressdT                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2IsochoricStressdFe,                      floatMatrix, setdPK2IsochoricStressdFe                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2StressdT,                                floatVector, setdPK2StressdT                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdT,                             floatVector, setdCauchyStressdT                             )

        };

    }

}

#endif
