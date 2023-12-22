/**
  ******************************************************************************
  * \file tardigrade_hydraThermalExpansion.h
  ******************************************************************************
  * An implementation of thermal expansion using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_THERMAL_EXPANSION_H
#define TARDIGRADE_HYDRA_THERMAL_EXPANSION_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace thermalExpansion{

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
         * A class which defines a thermal expansion residual
         */
        class residual : public tardigradeHydra::residualBase{

            public:

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int thermalConfigurationIndex, const floatVector &parameters ) : tardigradeHydra::residualBase( hydra, numEquations ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &thermalConfigurationIndex: The index of the thermal
                     *     configuration.
                     * \param &parameters: The parameters for the model as required by
                     *     constitutiveTools::quadraticThermalExpansion
                     */

                    _thermalConfigurationIndex = thermalConfigurationIndex;

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameters( parameters ) );

                }

                void setReferenceTemperature( const floatType &referenceTemperature );

                void setLinearParameters( const floatVector &linearParameters );

                void setQuadraticParameters( const floatVector &quadraticParameters );

                //! Get the index of the thermal deformation configuration
                const unsigned int* getThermalConfigurationIndex( ){ return &_thermalConfigurationIndex; };

                //! Get the reference temperature for the thermal expansion
                const floatType* getReferenceTemperature( ){ return &_referenceTemperature; }

                //! Get the linear parameters for the thermal expansion
                const floatVector* getLinearParameters( ){ return &_linearParameters; }

                //! Get the quadratic parameters for the thermal expansion
                const floatVector* getQuadraticParameters( ){ return &_quadraticParameters; }

            protected:

                virtual void setThermalGreenLagrangeStrain( );

                virtual void setThermalDeformationGradient( );

                virtual void setdThermalGreenLagrangeStraindT( );

                virtual void setdThermalDeformationGradientdT( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

            private:

                // Friend classes
                friend class tardigradeHydra::thermalExpansion::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBase::residualBase;

                using tardigradeHydra::residualBase::setResidual;

                using tardigradeHydra::residualBase::setJacobian;

                using tardigradeHydra::residualBase::setdRdF;

                using tardigradeHydra::residualBase::setdRdT;

                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                unsigned int _thermalConfigurationIndex;

                floatType _referenceTemperature;

                floatVector _linearParameters;

                floatVector _quadraticParameters;

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, thermalGreenLagrangeStrain,    floatVector, setThermalGreenLagrangeStrain    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dThermalGreenLagrangeStraindT, floatVector, setdThermalGreenLagrangeStraindT )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, thermalDeformationGradient,    floatVector, setThermalDeformationGradient    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dThermalDeformationGradientdT, floatVector, setdThermalDeformationGradientdT )

                virtual void decomposeParameters( const floatVector &parameters );

        };

    }

}

#endif
