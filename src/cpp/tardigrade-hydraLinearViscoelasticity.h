/**
  ******************************************************************************
  * \file tardigrade-hydraLinearViscoelasticity.h
  ******************************************************************************
  * An implementation of linear viscoelasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_LINEAR_VISCOELASTICITY_H
#define TARDIGRADE_HYDRA_LINEAR_VISCOELASTICITY_H

#define USE_EIGEN
#include<vector_tools.h>
#include<tardigrade-hydraLinearElasticity.h>

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
        const std::string __BASENAME__ = file_name(__FILE__);
        const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of("."));
    
        typedef errorTools::Node errorNode; //!< Redefinition for the error node
        typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
        typedef double floatType; //!< Define the float values type.
        typedef std::vector< floatType > floatVector; //!< Define a vector of floats
        typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats 
    
        class residual : public tardigradeHydra::linearElasticity::residual {
            /*!
             * A residual class for a linear-elastic material model where the stress is computed
             * in the reference configuration and pushed forward to the current configuration.
             */
        
            public:

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const floatVector &parameters,
                          const unsigned int viscoelasticISVLowerIndex,
                          const unsigned int viscoelasticISVUpperIndex,
                          const floatType integrationAlpha=0. ) : tardigradeHydra::linearElasticity::residual( hydra, numEquations ), _viscoelasticISVLowerIndex( viscoelasticISVLowerIndex ), _viscoelasticISVUpperIndex( viscoelasticISVUpperIndex ), _integrationAlpha( integrationAlpha ){
    
                    ERROR_TOOLS_CATCH( decomposeParameterVector( parameters ) );
    
                }

                const floatVector* getUpdatedStateVariables( ); //TODO: This will eventually need to be something in the base class

                //! Get the lower index of the viscoelastic ISVs from the non-nonlinear solve state variable vector
                const unsigned int* getViscoelasticISVLowerIndex( ){ return &_viscoelasticISVLowerIndex; }

                //! Get the upper (but not including) index of the viscoelastic ISVs from the non-nonlinear solve state variable vector
                const unsigned int* getViscoelasticISVUpperIndex( ){ return &_viscoelasticISVUpperIndex; }

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

                virtual void decomposeDeformation( const floatVector &Fe, floatType &Je, floatVector &Fehat );

                const floatVector* getVolumetricTemperatureParameters( ){ return &_volumetricTemperatureParameters; }

                const floatVector* getIsochoricTemperatureParameters( ){ return &_isochoricTemperatureParameters; }

                virtual void decomposeElasticDeformation( );

                virtual void decomposePreviousElasticDeformation( );

                void setJe( const floatType &Je );

                const floatType* getJe( );

                void setFehat( const floatVector &Fehat );

                const floatVector* getFehat( );

                void setPreviousJe( const floatType &Je );

                const floatType* getPreviousJe( );

                void setPreviousFehat( const floatVector &Fehat );

                const floatVector* getPreviousFehat( );

                virtual void setdJedFe( );

                void setdJedFe( const floatVector &dJedFe );

                const floatVector *getdJedFe( );

                virtual void setdFehatdFe( );

                void setdFehatdFe( const floatMatrix &dFehatdFe );

                const floatMatrix *getdFehatdFe( );

                void setNumStateVariables( const unsigned int numStateVariables );

                unsigned int *getNumStateVariables( ){ return &_numStateVariables; }

                virtual void decomposeStateVariableVector( floatVector &volumetricStateVariables, floatVector &isochoricStateVariables );

                virtual floatType computeRateMultiplier( const floatVector &variables, const floatVector &parameters );

                virtual void setVolumetricRateMultiplier( const floatType &rateMultiplier );

                virtual void setPreviousVolumetricRateMultiplier( const floatType &previousRateMultiplier );

                const floatType* getVolumetricRateMultiplier( );

                const floatType* getPreviousVolumetricRateMultiplier( );

                virtual void setIsochoricRateMultiplier( const floatType &rateMultiplier );

                virtual void setPreviousIsochoricRateMultiplier( const floatType &previousRateMultiplier );

                const floatType* getIsochoricRateMultiplier( );

                const floatType* getPreviousIsochoricRateMultiplier( );

                void setVolumetricTemperatureParameters( const floatVector &parameters );

                void setIsochoricTemperatureParameters( const floatVector &parameters );

            protected:

                virtual void setNumVolumetricViscousTerms( const unsigned int &num );

                virtual void setNumIsochoricViscousTerms( const unsigned int &num );

                virtual void setKinf( const floatType &Kinf );

                virtual void setGinf( const floatType &Ginf );

                virtual void setVolumetricModuli( const floatVector &Ks );

                virtual void setIsochoricModuli( const floatVector &Gs );
 
                virtual void setVolumetricTaus( const floatVector &taus );

                virtual void setIsochoricTaus( const floatVector &taus );
 
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
        
                dataStorage< floatType > _Je;

                dataStorage< floatVector > _Fehat;

                dataStorage< floatType > _previousJe;

                dataStorage< floatVector > _previousFehat;

                dataStorage< floatVector > _dJedFe;

                dataStorage< floatMatrix > _dFehatdFe;

                dataStorage< floatType > _volumetricRateMultiplier;

                dataStorage< floatType > _previousVolumetricRateMultiplier;

                dataStorage< floatType > _isochoricRateMultiplier;

                dataStorage< floatType > _previousIsochoricRateMultiplier;

                dataStorage< floatVector > _currentStateVariables;

                virtual void setPK2Stress( ) override;
    
                virtual void decomposeParameterVector( const floatVector &parameters );

                void setVolumetricRateMultiplier( );

                void setPreviousVolumetricRateMultiplier( );

                void setIsochoricRateMultiplier( );

                void setPreviousIsochoricRateMultiplier( );
        };

    }

}

#endif
