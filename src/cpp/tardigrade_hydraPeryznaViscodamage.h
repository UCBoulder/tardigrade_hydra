/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaViscodamage.h
  ******************************************************************************
  * An implementation of peryznaViscodamage using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_PERYZNA_VISCODAMAGE_H
#define TARDIGRADE_HYDRA_PERYZNA_VISCODAMAGE_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydraPeryznaViscoplasticity.h>

namespace tardigradeHydra{

    namespace peryznaViscodamage{

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
         * A class which defines a viscodamage residual
         */
        class residual : public tardigradeHydra::peryznaViscoplasticity::residual{

            public:

                //! Default constructor
                residual( ) : tardigradeHydra::peryznaViscoplasticity::residual( ), _elasticConfigurationIndex( 0 ){ }

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const unsigned int &elasticConfigurationIndex = 0, const floatType integrationParameter = 0.5 ) : tardigradeHydra::peryznaViscoplasticity::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters, integrationParameter ),
                    _elasticConfigurationIndex( elasticConfigurationIndex ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &damageConfigurationIndex: The index of the damage
                     *     configuration.
                     * \param &stateVariableIndices: The indices of the plastic state variables
                     * \param &parameters: The parameters for the model
                     * \param &elasticConfigurationIndex: The index of the elastic configuration
                     * \param &integrationParameter: The integration parameter for the function. 0 is explicit, 1 is implicit.
                     */

                    //Setting the contribution of damage to the calculation of the drag stress to zero
                    set_dragStressParameters( { parameters[ 1 ], parameters[  2 ], 0. } );

                    //Setting the contribution of damage to the hardening of the damage state variable to zero
                    set_hardeningParameters(  { parameters[ 9 ], parameters[ 10 ], 0. } );

                }

            public:

                 // Friend classes
                friend class tardigradeHydra::peryznaViscoplasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::peryznaViscoplasticity::residual::residual;

                using tardigradeHydra::peryznaViscoplasticity::residual::setResidual;

                using tardigradeHydra::peryznaViscoplasticity::residual::setJacobian;

                using tardigradeHydra::peryznaViscoplasticity::residual::setdRdF;

                using tardigradeHydra::peryznaViscoplasticity::residual::setdRdT;

                using tardigradeHydra::peryznaViscoplasticity::residual::setAdditionalDerivatives;

                //! Get the configuration of the damage
                const unsigned int *getDamageConfigurationIndex( ){ return &( *getPlasticConfigurationIndex( ) ); }

                virtual void setDamage( );

                virtual void setDamageJacobians( );

                virtual void setAllDamageJacobians( );

                virtual void setDamageJacobians( const bool withPrevious );

                virtual void setDamageDeformationGradient( );

                virtual void setDamageDeformationGradientJacobians( );

                virtual void setAllDamageDeformationGradientJacobians( );

                virtual void setDamageDeformationGradientJacobians( const bool withPrevious );

                //! Get the index of the elastic configuration
                const unsigned int *getElasticConfigurationIndex( ){ return &_elasticConfigurationIndex; }

            private:

                unsigned int _elasticConfigurationIndex;

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, damage,                                            floatType,   setDamage                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, damageDeformationGradient,                         floatVector, setDamageDeformationGradient             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedCauchyStress,                              floatVector, setDamageJacobians                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedF,                                         floatVector, setDamageJacobians                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedSubFs,                                     floatVector, setDamageJacobians                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedT,                                         floatType,   setDamageJacobians                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedStateVariables,                            floatVector, setDamageJacobians                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedPreviousCauchyStress,                      floatVector, setAllDamageJacobians                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedPreviousF,                                 floatVector, setAllDamageJacobians                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedPreviousSubFs,                             floatVector, setAllDamageJacobians                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedPreviousT,                                 floatType,   setAllDamageJacobians                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamagedPreviousStateVariables,                    floatVector, setAllDamageJacobians                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdCauchyStress,           floatMatrix, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdF,                      floatMatrix, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdSubFs,                  floatMatrix, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdT,                      floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdStateVariables,         floatMatrix, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousCauchyStress,   floatMatrix, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousF,              floatMatrix, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousSubFs,          floatMatrix, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousT,              floatVector, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousStateVariables, floatMatrix, setAllDamageDeformationGradientJacobians )

        };

    }

}

#endif
