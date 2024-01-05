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

    namespace micromorphicDruckerPragerPlasticity{

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

        /*!
         * The residual for a micromorphic Drucker Prager plasticity model
         */
        class residual : public tardigradeHydra::residualBaseMicromorphic {

            public:

                using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

                residual( hydraBaseMicromorphic *_hydra, const unsigned int &_numEquations, const unsigned int &plasticConfigurationIndex,
                          const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const floatType integrationParameter = 0.5 )
                        : tardigradeHydra::residualBaseMicromorphic( _hydra, _numEquations ){
                    /*!
                     * The main initialization constructor for the Drucker Prager plasticity residual
                     * 
                     * \param *_hydra: A pointer to the containing hydra class
                     * \param &_numEquations: The number of equations the residual defines
                     * \param &plasticConfigurationIndex: The index of the configuration which represents the plastic deformation
                     * \param &stateVariableIndices: The indices of the plastic state variables
                     * \param &parameters: The parameter vector
                     * \param &integrationParameter: The integration parameter for the function. 0 is explicit, 1 is implicit.
                     */

                    _plasticConfigurationIndex = plasticConfigurationIndex;

                    _stateVariableIndices = stateVariableIndices;

                    _integrationParameter = integrationParameter;

                }

                const unsigned int* getPlasticConfigurationIndex( );

                const std::vector< unsigned int >* getStateVariableIndices( );

                const floatType* getIntegrationParameter( );

            protected:

                virtual void setMacroDrivingStress( );

                virtual void setSymmetricMicroDrivingStress( );

                virtual void setHigherOrderDrivingStress( );

                virtual void setPreviousMacroDrivingStress( );

                virtual void setPreviousSymmetricMicroDrivingStress( );

                virtual void setPreviousHigherOrderDrivingStress( );

                virtual void setDrivingStresses( const bool isPrevious );

                virtual void setdMacroDrivingStressdMacroStress( );

                virtual void setdSymmetricMicroDrivingStressdMicroStress( );

                virtual void setdHigherOrderDrivingStressdHigherOrderStress( );

                virtual void setdMacroDrivingStressdF( );

                virtual void setdSymmetricMicroDrivingStressdF( );

                virtual void setdHigherOrderDrivingStressdF( );

                virtual void setdHigherOrderDrivingStressdChi( );

                virtual void setdMacroDrivingStressdFn( );

                virtual void setdSymmetricMicroDrivingStressdFn( );

                virtual void setdHigherOrderDrivingStressdFn( );

                virtual void setdHigherOrderDrivingStressdChin( );

                virtual void setDrivingStressesJacobians( const bool isPrevious );

            private:

                unsigned int _plasticConfigurationIndex; //! The index of the plastic configuration

                std::vector< unsigned int > _stateVariableIndices; //! The indices of the state variables in the global solve

                floatType _integrationParameter; //! The integration parameter (0 is explicit, 1 is implicit)

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroDrivingStress,                          floatVector, setMacroDrivingStress                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, symmetricMicroDrivingStress,                 floatVector, setSymmetricMicroDrivingStress                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, higherOrderDrivingStress,                    floatVector, setHigherOrderDrivingStress                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdMacroStress,             floatMatrix, setdMacroDrivingStressdMacroStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdMicroStress,    floatMatrix, setdSymmetricMicroDrivingStressdMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdHigherOrderStress, floatMatrix, setdHigherOrderDrivingStressdHigherOrderStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdF,                       floatMatrix, setdMacroDrivingStressdF                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdF,              floatMatrix, setdSymmetricMicroDrivingStressdF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdF,                 floatMatrix, setdHigherOrderDrivingStressdF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChi,               floatMatrix, setdHigherOrderDrivingStressdChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdFn,                      floatMatrix, setdMacroDrivingStressdFn                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdFn,             floatMatrix, setdSymmetricMicroDrivingStressdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdFn,                floatMatrix, setdHigherOrderDrivingStressdFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChin,              floatMatrix, setdHigherOrderDrivingStressdChin              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroDrivingStress,                  floatVector, setPreviousMacroDrivingStress                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSymmetricMicroDrivingStress,         floatVector, setPreviousSymmetricMicroDrivingStress         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHigherOrderDrivingStress,            floatVector, setPreviousHigherOrderDrivingStress            )

        };

    }

}

#endif
