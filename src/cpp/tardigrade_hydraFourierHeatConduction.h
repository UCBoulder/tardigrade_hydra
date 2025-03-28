/**
  ******************************************************************************
  * \file tardigrade_hydraFourierHeatConduction.h
  ******************************************************************************
  * An implementation of Fourier heat conduction using the hydra framework.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_FOURIER_HEAT_CONDUCTION_H
#define TARDIGRADE_HYDRA_FOURIER_HEAT_CONDUCTION_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace fourierHeatConduction{

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

        /*!
         * A residual class for a linear-elastic material model where the stress is computed
         * in the reference configuration and pushed forward to the current configuration.
         */
        class residual : public tardigradeHydra::residualBase {
        
            public:

                /*!
                 * The main initialization constructor for the linear elastic residual
                 * 
                 * \param *hydra: A pointer to the containing hydra class
                 * \param &numEquations: The number of equations the residual defines
                 * \param &parameters: The parameter vector
                 * \param &temperatureGradientIndex: The index in the additional DOF array where the temperature gradient exists
                 * \param &expectedParameterVectorSize: The expected size of the parameter vector (defaults to 1)
                 * \param &heatFluxIndex: The index of the heat flux in the residual vector (defaults to 17)
                 */
                residual(
                    tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const floatVector &parameters,
                    const unsigned int temperatureGradientIndex,
                    const unsigned int expectedParameterVectorSize = 1,
                    const unsigned int heatFluxIndex = 17
                ) : tardigradeHydra::residualBase( hydra, numEquations ),
                _expectedParameterVectorSize( expectedParameterVectorSize ),
                _temperatureGradientIndex( temperatureGradientIndex ),
                _heatFluxIndex( heatFluxIndex ){
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameterVector( parameters ) );

                }

            //! Get the expected parameter vector size
            auto get_expectedParameterVectorSize( ){ return _expectedParameterVectorSize; }

            //! Get the index in the additional DOF vector where the temperature gradient is stored
            auto get_temperatureGradientIndex( ){ return _temperatureGradientIndex; }

            //! Get the index in the residual vector where the heat flux is stored
            auto get_heatFluxIndex( ){ return _heatFluxIndex; }

            const floatVector *get_conductivityParameters( ){
                /*!
                 * Get the conductivity parameters
                 */

                return &_conductivityParameters;

            }

            protected:

                virtual void setConductivityParameters(
                    const floatVector &conductivityParameters
                ){
                    /*!
                     * Set the values of the conductivity parameters
                     * 
                     * \param &conductivityParameters: The values of the conductivity parameters
                     */

                     _conductivityParameters = conductivityParameters;

                }

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdAdditionalDOF( ) override;

                virtual void decomposeParameterVector( const floatVector &parameters );

                virtual void setTemperatureGradient( );

                virtual void setPreviousTemperatureGradient( );

                virtual void setConductivity( );

                virtual void setPreviousConductivity( );

                virtual void setHeatFlux( );

                virtual void setPreviousHeatFlux( );

                virtual void setdHeatFluxdGradT( );
     
            private:
    
                // Friend classes
                friend class tardigradeHydra::fourierHeatConduction::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!
        
                using tardigradeHydra::residualBase::residualBase;
        
                using tardigradeHydra::residualBase::setResidual;
        
                using tardigradeHydra::residualBase::setJacobian;
        
                using tardigradeHydra::residualBase::setdRdF;
        
                using tardigradeHydra::residualBase::setdRdT;
        
                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, temperatureGradient,               floatVector,          setTemperatureGradient )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, previousTemperatureGradient,       floatVector,  setPreviousTemperatureGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, conductivity,                        floatType,                 setConductivity )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousConductivity,                floatType,         setPreviousConductivity )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, heatFlux,                          floatVector,                     setHeatFlux )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHeatFlux,                  floatVector,             setPreviousHeatFlux )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHeatFluxdGradT,             secondOrderTensor,              setdHeatFluxdGradT )

                floatVector _conductivityParameters; //!< The parameters used in the calculation of the conductivity

                const unsigned int _expectedParameterVectorSize; //!< The expected size of the parameter vector

                const unsigned int _temperatureGradientIndex; //!< The index in the additional DOF vector where the temperature gradient is located

                const unsigned int _heatFluxIndex; //!< The index of the heat flux in the output residual vector

        };

    }

}

#endif
