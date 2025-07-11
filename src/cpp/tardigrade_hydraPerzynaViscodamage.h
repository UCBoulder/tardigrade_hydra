/**
  ******************************************************************************
  * \file tardigrade_hydraPerzynaViscodamage.h
  ******************************************************************************
  * An implementation of perzynaViscodamage using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_PERZYNA_VISCODAMAGE_H
#define TARDIGRADE_HYDRA_PERZYNA_VISCODAMAGE_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydraPerzynaViscoplasticity.h>

namespace tardigradeHydra{

    namespace perzynaViscodamage{

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
        class residual : public tardigradeHydra::perzynaViscoplasticity::residual{

            public:

                //! Default constructor
                residual( ) : tardigradeHydra::perzynaViscoplasticity::residual( ), _elasticConfigurationIndex( 0 ){ }

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const unsigned int &elasticConfigurationIndex = 0, const floatType integrationParameter = 0.5 ) : tardigradeHydra::perzynaViscoplasticity::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters, integrationParameter ),
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
                friend class tardigradeHydra::perzynaViscoplasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::perzynaViscoplasticity::residual::residual;

                using tardigradeHydra::perzynaViscoplasticity::residual::setResidual;

                using tardigradeHydra::perzynaViscoplasticity::residual::setJacobian;

                using tardigradeHydra::perzynaViscoplasticity::residual::setdRdF;

                using tardigradeHydra::perzynaViscoplasticity::residual::setdRdT;

                using tardigradeHydra::perzynaViscoplasticity::residual::setAdditionalDerivatives;

                using tardigradeHydra::perzynaViscoplasticity::residual::setStateVariableEvolutionRates;

                //! Get the configuration of the damage
                const unsigned int *getDamageConfigurationIndex( ){ return &( *getPlasticConfigurationIndex( ) ); }

                virtual void setStateVariableEvolutionRates( const bool isPrevious ) override;

                virtual void setStateVariableEvolutionRateDerivatives( const bool isPrevious ) override;

                virtual void setDamage( );

                virtual void setDamageJacobians( );

                virtual void setAllDamageJacobians( );

                virtual void setDamageJacobians( const bool withPrevious );

                virtual void setDamageDeformationGradient( );

                virtual void setDamageDeformationGradientJacobians( );

                virtual void setAllDamageDeformationGradientJacobians( );

                virtual void setDamageDeformationGradientJacobians( const bool withPrevious );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdT( ) override;

                //! Get the index of the elastic configuration
                const unsigned int *getElasticConfigurationIndex( ){ return &_elasticConfigurationIndex; }

                //! The index of the scalar damage
                unsigned int damageISVIndex = 1;

                virtual void addParameterizationInfo( std::string &parameterization_info ) override;

            protected:

                virtual void decomposeParameters( const floatVector &parameters ) override;

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

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdCauchyStress,           floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdF,                      floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdSubFs,                  floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdT,                      floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdStateVariables,         floatVector, setDamageDeformationGradientJacobians    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousCauchyStress,   floatVector, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousF,              floatVector, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousSubFs,          floatVector, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousT,              floatVector, setAllDamageDeformationGradientJacobians )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDamageDeformationGradientdPreviousStateVariables, floatVector, setAllDamageDeformationGradientJacobians )

        };

    }

}

#endif
