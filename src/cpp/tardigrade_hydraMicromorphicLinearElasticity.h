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
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace linearElasticity{

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

    }

}

#endif
