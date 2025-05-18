/**
  ******************************************************************************
  * \file tardigrade_hydraLinearInternalEnergy.h
  ******************************************************************************
  * An implementation of Fourier heat conduction using the hydra framework.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_LINEAR_INTERNAL_ENERGY_H
#define TARDIGRADE_HYDRA_LINEAR_INTERNAL_ENERGY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace linearInternalEnergy{

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
                 * The main initialization constructor for the linear temperature-internal energy relation
                 * 
                 * \param *hydra: A pointer to the containing hydra class
                 * \param &numEquations: The number of equations the residual defines
                 * \param &parameters: The parameter vector
                 * \param &expectedParameterVectorSize: The expected size of the parameter vector (defaults to 1)
                 * \param &internalEnergyIndex: The index of the predicted internal energy residual vector (defaults to 9)
                 */
                residual(
                    tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const floatVector &parameters,
                    const unsigned int expectedParameterVectorSize = 1,
                    const unsigned int internalEnergyIndex = 9
                ) : tardigradeHydra::residualBase( hydra, numEquations ),
                _expectedParameterVectorSize( expectedParameterVectorSize ),
                _internalEnergyIndex( internalEnergyIndex ){
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameterVector( parameters ) );

                }

            //! Get the expected parameter vector size
            auto get_expectedParameterVectorSize( ){ return _expectedParameterVectorSize; }

            //! Get the index in the residual vector where the internal energy is stored
            auto get_internalEnergyIndex( ){ return _internalEnergyIndex; }

            const floatVector *get_specificHeatParameters( ){
                /*!
                 * Get the specific heat parameters
                 */

                return &_specificHeatParameters;

            }

            protected:

                virtual void setSpecificHeatParameters(
                    const floatVector &specificHeatParameters
                ){
                    /*!
                     * Set the values of the specific heat parameters
                     * 
                     * \param &specificHeatParameters: The values of the specific heat parameters
                     */

                     _specificHeatParameters = specificHeatParameters;

                }

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdAdditionalDOF( ) override;

                virtual void decomposeParameterVector( const floatVector &parameters );

                virtual void setSpecificHeat( );

                virtual void setPreviousSpecificHeat( );

                virtual void setInternalEnergy( );

                virtual void setPreviousInternalEnergy( );

                virtual void setdInternalEnergydT( );
     
            private:
    
                // Friend classes
                friend class tardigradeHydra::linearInternalEnergy::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!
        
                using tardigradeHydra::residualBase::residualBase;
        
                using tardigradeHydra::residualBase::setResidual;
        
                using tardigradeHydra::residualBase::setJacobian;
        
                using tardigradeHydra::residualBase::setdRdF;
        
                using tardigradeHydra::residualBase::setdRdT;
        
                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, specificHeat,           floatType,           setSpecificHeat )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSpecificHeat,   floatType,   setPreviousSpecificHeat )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, internalEnergy,         floatType,         setInternalEnergy )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, previousInternalEnergy, floatType, setPreviousInternalEnergy )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dInternalEnergydT,      floatType,      setdInternalEnergydT )

                floatVector _specificHeatParameters; //!< The parameters used in the calculation of the specific heat

                unsigned int _expectedParameterVectorSize; //!< The expected size of the parameter vector

                unsigned int _internalEnergyIndex; //!< The index of the internal energy in the output residual vector

        };

    }

}

#endif
