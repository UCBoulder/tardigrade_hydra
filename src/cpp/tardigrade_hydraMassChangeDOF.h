/**
  ******************************************************************************
  * \file tardigrade_hydraMassChangeDOF.h
  ******************************************************************************
  * An implementation of the mass-change residual where the mass change velocity
  * gradient is defined via a DOF.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_MASS_CHANGE_DOF_H
#define TARDIGRADE_HYDRA_MASS_CHANGE_DOF_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace massChangeDOF{

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
         * A class which defines a mass-change residual where the mass change velocity gradient
         * is defined by a degree of freedom gradient located in the additional DOF vector
         */
        class residual : public tardigradeHydra::residualBase{

            public:

                residual(
                    tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                    const unsigned int massChangeConfigurationIndex, const unsigned int massChangeVelocityGradientIndex,
                    const floatVector &parameters, const floatType integrationParameter = 0.5
                ) : tardigradeHydra::residualBase( hydra, numEquations ), _integrationParameter( integrationParameter ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &massChangeConfigurationIndex: The index of the mass-change configuration
                     * \param &massChangeVelocityGradientIndex: The index of the current configuration mass-change velocity gradient in the additional dof vector
                     * \param &parameters: The parameters for the model
                     * \param integrationParameter: The parameter of the integration 0 is explicit, 1 is implicit
                     */

                    _massChangeConfigurationIndex = massChangeConfigurationIndex;

                    _massChangeVelocityGradientIndex = massChangeVelocityGradientIndex;

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeAdditionalDOF( ) );

                }

                //! Get the index of the mass-change configuration
                const unsigned int getMassChangeConfigurationIndex( ){ return _massChangeConfigurationIndex; }

                //! Get the index of the velocity gradient in the DOF vector
                const unsigned int getMassChangeVelocityGradientIndex( ){ return _massChangeVelocityGradientIndex; }

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

                virtual void setMassChangeIntermediateVelocityGradient( const bool &isPrevious );

                virtual void setMassChangeIntermediateVelocityGradientDerivatives( const bool &isPrevious );

                virtual void setMassChangeIntermediateVelocityGradient( );

                virtual void setdMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient( );

                virtual void setdMassChangeIntermediateVelocityGradientdDeformationGradient( );

                virtual void setdMassChangeIntermediateVelocityGradientdSubDeformationGradients( );

                virtual void setPreviousMassChangeIntermediateVelocityGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( );

                virtual void setMassChangeDeformationGradient( );

                virtual void setMassChangeDeformationGradientDerivatives( const bool &computePrevious );

                virtual void setdMassChangeDeformationGradientdMassChangeVelocityGradient( );

                virtual void setdMassChangeDeformationGradientdDeformationGradient( );

                virtual void setdMassChangeDeformationGradientdSubDeformationGradients( );

                virtual void setdMassChangeDeformationGradientdPreviousMassChangeVelocityGradient( );

                virtual void setdMassChangeDeformationGradientdPreviousDeformationGradient( );

                virtual void setdMassChangeDeformationGradientdPreviousSubDeformationGradients( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdAdditionalDOF( ) override;

            private:

                // Friend classes
                friend class tardigradeHydra::massChangeDOF::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBase::residualBase;

                using tardigradeHydra::residualBase::setResidual;

                using tardigradeHydra::residualBase::setJacobian;

                using tardigradeHydra::residualBase::setdRdF;

                using tardigradeHydra::residualBase::setdRdT;

                using tardigradeHydra::residualBase::setdRdAdditionalDOF;

                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                unsigned int _massChangeConfigurationIndex;

                unsigned int _massChangeVelocityGradientIndex;

                floatType _integrationParameter;

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(
                    private, massChangeVelocityGradient,
                    dimVector, unexpectedError
                )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(
                    private, previousMassChangeVelocityGradient,
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
                    private,              massChangeIntermediateVelocityGradient,
                    secondOrderTensor, setMassChangeIntermediateVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient,
                    fourthOrderTensor, setdMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dMassChangeIntermediateVelocityGradientdDeformationGradient,
                    fourthOrderTensor, setdMassChangeIntermediateVelocityGradientdDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dMassChangeIntermediateVelocityGradientdSubDeformationGradients,
                    floatVector, setdMassChangeIntermediateVelocityGradientdSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              previousMassChangeIntermediateVelocityGradient,
                    secondOrderTensor, setPreviousMassChangeIntermediateVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient,
                    fourthOrderTensor, setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient,
                    fourthOrderTensor, setdPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients,
                    floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              massChangeDeformationGradient,
                    secondOrderTensor, setMassChangeDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dMassChangeDeformationGradientdMassChangeVelocityGradient,
                    fourthOrderTensor, setdMassChangeDeformationGradientdMassChangeVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,              dMassChangeDeformationGradientdDeformationGradient,
                    fourthOrderTensor, setdMassChangeDeformationGradientdDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(
                    private,        dMassChangeDeformationGradientdSubDeformationGradients,
                    floatVector, setdMassChangeDeformationGradientdSubDeformationGradients
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dMassChangeDeformationGradientdPreviousMassChangeVelocityGradient,
                    fourthOrderTensor, setdMassChangeDeformationGradientdPreviousMassChangeVelocityGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,              dMassChangeDeformationGradientdPreviousDeformationGradient,
                    fourthOrderTensor, setdMassChangeDeformationGradientdPreviousDeformationGradient
                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(
                    private,        dMassChangeDeformationGradientdPreviousSubDeformationGradients,
                    floatVector, setdMassChangeDeformationGradientdPreviousSubDeformationGradients
                )

        };

    }

}

#endif
