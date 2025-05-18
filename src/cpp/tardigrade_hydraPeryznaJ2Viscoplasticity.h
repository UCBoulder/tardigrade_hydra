/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaJ2Viscoplasticity.h
  ******************************************************************************
  * An implementation of peryznaJ2Viscoplasticity using the hydra framework.
  * Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_PERYZNA_J2_VISCOPLASTICITY_H
#define TARDIGRADE_HYDRA_PERYZNA_J2_VISCOPLASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydraPeryznaViscoplasticity.h>

namespace tardigradeHydra{

    namespace peryznaJ2Viscoplasticity{

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

        template <typename T> int sgn(T val){
            /*!
             * Get the sign of the value
             * 
             * \param val: The value to compute the sign of
             */
            return (T(0) < val) - (val < T(0));
        }

        /*!
         * A class which defines a J2 viscoplastic residual
         */
        class residual : public tardigradeHydra::peryznaViscoplasticity::residual{

            public:

                //! Default constructor
                residual( ) : tardigradeHydra::peryznaViscoplasticity::residual( ), _elasticConfigurationIndex( 0 ){ }

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const unsigned int &elasticConfigurationIndex = 0, const floatType integrationParameter = 0.5 ) : tardigradeHydra::peryznaViscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters, integrationParameter ),
                    _elasticConfigurationIndex( elasticConfigurationIndex ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &plasticConfigurationIndex: The index of the plastic
                     *     configuration.
                     * \param &stateVariableIndices: The indices of the plastic state variables
                     * \param &parameters: The parameters for the model
                     * \param &elasticConfigurationIndex: The index of the elastic configuration
                     * \param &integrationParameter: The integration parameter for the function. 0 is explicit, 1 is implicit.
                     */

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

                //! Get the index of the elastic configuration
                const unsigned int *getElasticConfigurationIndex( ){ return &_elasticConfigurationIndex; }

                //! The index of the scalar damage
                unsigned int damageISVIndex = 1;

            protected:

                using tardigradeHydra::peryznaViscoplasticity::residual::setHardeningFunction;

                using tardigradeHydra::peryznaViscoplasticity::residual::setDragStress;

                using tardigradeHydra::peryznaViscoplasticity::residual::setYieldFunction;

                virtual void setYieldFunction( const bool isPrevious ) override;

                virtual void setYieldFunctionDerivatives( const bool isPrevious ) override;

                virtual void setSignTerm( );

                virtual void setPreviousSignTerm( );

                virtual void setDragStress( const bool isPrevious ) override;

                virtual void setDragStressDerivatives( const bool isPrevious ) override;

                virtual void setHardeningFunction( const bool isPrevious ) override;

                virtual void setHardeningFunctionDerivatives( const bool isPrevious ) override;

                virtual void decomposeParameters( const floatVector &parameters ) override;

            private:

                unsigned int _elasticConfigurationIndex;

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, signTerm,         floatType,   setSignTerm         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSignTerm, floatType,   setPreviousSignTerm )


        };

    }

}

#endif
