/**
  ******************************************************************************
  * \file tardigrade_hydraDOFVelocityGradientDeformation.h
  ******************************************************************************
  * An implementation of the evolution of deformation where the velocity
  * gradient is defined via a DOF.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_DOF_VELOCITY_GRADIENT_DEFORMATION_H
#define TARDIGRADE_HYDRA_DOF_VELOCITY_GRADIENT_DEFORMATION_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace dofVelocityGradientDeformation{

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
         * A class which defines a deformation where the velocity gradient of that deformation is
         * defined by the additional dof vector
         */
        class residual : public tardigradeHydra::residualBase{

            public:

                residual(
                    tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                    const unsigned int dofConfigurationIndex, const unsigned int dofVelocityGradientIndex,
                    const floatVector &parameters, const floatType integrationParameter = 0.5
                ) : tardigradeHydra::residualBase( hydra, numEquations ), _integrationParameter( integrationParameter ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &dofConfigurationIndex: The index of the mass-change configuration
                     * \param &dofVelocityGradientIndex: The index of the current configuration velocity gradient in the additional dof vector
                     * \param &parameters: The parameters for the model
                     * \param integrationParameter: The parameter of the integration 0 is explicit, 1 is implicit
                     */

                    _dofConfigurationIndex = dofConfigurationIndex;

                    _dofVelocityGradientIndex = dofVelocityGradientIndex;

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeAdditionalDOF( ) );

                }

                //! Get the index of the mass-change configuration
                const unsigned int getDOFConfigurationIndex( ){ return _dofConfigurationIndex; }

                //! Get the index of the velocity gradient in the DOF vector
                const unsigned int getDOFVelocityGradientIndex( ){ return _dofVelocityGradientIndex; }

                //! Get the integration parameter 0 for explicit, 1 for implicit
                const floatType getIntegrationParameter( ){ return _integrationParameter; }

                virtual void suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                          std::vector< floatType > &values ) override;

            protected:

                virtual void decomposeAdditionalDOF( );

                virtual void setPrecedingDeformationGradient( const bool &isPrevious );

                virtual void setPrecedingDeformationGradientDerivatives( const bool &isPrevious );

                virtual void setPrecedingDeformationGradient( );

                virtual void setPreviousPrecedingDeformationGradient( );

                virtual void setdPrecedingDeformationGradientdDeformationGradient( );

                virtual void setdPrecedingDeformationGradientdSubDeformationGradients( );

                virtual void setdPreviousPrecedingDeformationGradientdPreviousDeformationGradient( );

                virtual void setdPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( );

                virtual void setDOFIntermediateVelocityGradient( const bool &isPrevious );

                virtual void setDOFIntermediateVelocityGradientDerivatives( const bool &isPrevious );

                virtual void setDOFIntermediateVelocityGradient( );

                virtual void setdDOFIntermediateVelocityGradientdDOFVelocityGradient( );

                virtual void setdDOFIntermediateVelocityGradientdDeformationGradient( );

                virtual void setdDOFIntermediateVelocityGradientdSubDeformationGradients( );

                virtual void setPreviousDOFIntermediateVelocityGradient( );

                virtual void setdPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient( );

                virtual void setdPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient( );

                virtual void setdPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients( );

                virtual void setDOFDeformationGradient( );

                virtual void setDOFDeformationGradientDerivatives( const bool &computePrevious );

                virtual void setdDOFDeformationGradientdDOFVelocityGradient( );

                virtual void setdDOFDeformationGradientdDeformationGradient( );

                virtual void setdDOFDeformationGradientdSubDeformationGradients( );

                virtual void setdDOFDeformationGradientdPreviousDOFVelocityGradient( );

                virtual void setdDOFDeformationGradientdPreviousDeformationGradient( );

                virtual void setdDOFDeformationGradientdPreviousSubDeformationGradients( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdAdditionalDOF( ) override;

            private:

                // Friend classes
                friend class tardigradeHydra::dofVelocityGradientDeformation::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBase::residualBase;

                using tardigradeHydra::residualBase::setResidual;

                using tardigradeHydra::residualBase::setJacobian;

                using tardigradeHydra::residualBase::setdRdF;

                using tardigradeHydra::residualBase::setdRdT;

                using tardigradeHydra::residualBase::setdRdAdditionalDOF;

                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                unsigned int _dofConfigurationIndex;

                unsigned int _dofVelocityGradientIndex;

                floatType _integrationParameter;

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(
                    private, dofVelocityGradient,
                    dimVector, unexpectedError
                )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(
                    private, previousDOFVelocityGradient,
                    dimVector, unexpectedError
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              precedingDeformationGradient,
                    secondOrderTensor, setPrecedingDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dPrecedingDeformationGradientdDeformationGradient,
                    fourthOrderTensor, setdPrecedingDeformationGradientdDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dPrecedingDeformationGradientdSubDeformationGradients,
                    floatVector, setdPrecedingDeformationGradientdSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              previousPrecedingDeformationGradient,
                    secondOrderTensor, setPreviousPrecedingDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dPreviousPrecedingDeformationGradientdPreviousDeformationGradient,
                    fourthOrderTensor, setdPreviousPrecedingDeformationGradientdPreviousDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,        dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients,
                    floatVector, setdPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dofIntermediateVelocityGradient,
                    secondOrderTensor, setDOFIntermediateVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dDOFIntermediateVelocityGradientdDOFVelocityGradient,
                    fourthOrderTensor, setdDOFIntermediateVelocityGradientdDOFVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dDOFIntermediateVelocityGradientdDeformationGradient,
                    fourthOrderTensor, setdDOFIntermediateVelocityGradientdDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dDOFIntermediateVelocityGradientdSubDeformationGradients,
                    floatVector, setdDOFIntermediateVelocityGradientdSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              previousDOFIntermediateVelocityGradient,
                    secondOrderTensor, setPreviousDOFIntermediateVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient,
                    fourthOrderTensor, setdPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient,
                    fourthOrderTensor, setdPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients,
                    floatVector, setdPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dofDeformationGradient,
                    secondOrderTensor, setDOFDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dDOFDeformationGradientdDOFVelocityGradient,
                    fourthOrderTensor, setdDOFDeformationGradientdDOFVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dDOFDeformationGradientdDeformationGradient,
                    fourthOrderTensor, setdDOFDeformationGradientdDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dDOFDeformationGradientdSubDeformationGradients,
                    floatVector, setdDOFDeformationGradientdSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dDOFDeformationGradientdPreviousDOFVelocityGradient,
                    fourthOrderTensor, setdDOFDeformationGradientdPreviousDOFVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dDOFDeformationGradientdPreviousDeformationGradient,
                    fourthOrderTensor, setdDOFDeformationGradientdPreviousDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,        dDOFDeformationGradientdPreviousSubDeformationGradients,
                    floatVector, setdDOFDeformationGradientdPreviousSubDeformationGradients
                )

        };

    }

}

#endif
