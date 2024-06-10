/**
  ******************************************************************************
  * \file tardigrade_hydraMassChange.h
  ******************************************************************************
  * An implementation of the mass-change residual. Used as the basis for more
  * complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_MASS_CHANGE_BASE_H
#define TARDIGRADE_HYDRA_MASS_CHANGE_BASE_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace massChange{

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
         * A class which defines a mass-change residual
         *
         * Assumes that the additional DOF vector is comprised of
         * 
         * [density, mass_change_rate, mass_change_gradient_x, mass_change_gradient_y, mass_change_gradient_z]
         */
        class residual : public tardigradeHydra::residualBase{

            public:

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int massChangeConfigurationIndex, const floatVector &parameters, const floatType integrationParameter = 0.5 ) : tardigradeHydra::residualBase( hydra, numEquations ), _integrationParameter( integrationParameter ){
                    /*!
                     * The main constructor function
                     *
                     * \param *hydra: A reference to the containing hydra object
                     * \param &numEquations: The number of equations to be defined by
                     *     the residual
                     * \param &thermalConfigurationIndex: The index of the mass-change
                     * \param &parameters: The parameters for the model
                     */

                    _massChangeConfigurationIndex = massChangeConfigurationIndex;

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeAdditionalDOF( ) );

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameters( parameters ) );

                }

                const unsigned int *getMassChangeConfigurationIndex( ){ return &_massChangeConfigurationIndex; }

                const floatType *getIntegrationParameter( ){ return &_integrationParameter; }

                virtual void suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                          std::vector< floatType > &values ) override;

            protected:

                virtual void decomposeAdditionalDOF( );

                virtual void setMassChangeVelocityGradientTrace( const bool &isPrevious );

                virtual void setMassChangeVelocityGradientTraceDerivatives( const bool &isPrevious );

                virtual void setMassChangeVelocityGradientTrace( );

                virtual void setdMassChangeVelocityGradientTracedDensity( );

                virtual void setdMassChangeVelocityGradientTracedMassChangeRate( );

                virtual void setPreviousMassChangeVelocityGradientTrace( );

                virtual void setdPreviousMassChangeVelocityGradientTracedPreviousDensity( );

                virtual void setdPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( );

                virtual void setUnitDirectionVector( const bool &isPrevious );

                virtual void setUnitDirectionVectorDerivatives( const bool &isPrevious );

                virtual void setUnitDirectionVector( );

                virtual void setdUnitDirectionVectordMassChangeRateGradient( );

                virtual void setPreviousUnitDirectionVector( );

                virtual void setdPreviousUnitDirectionVectordPreviousMassChangeRateGradient( );

                virtual void setMassChangeVelocityGradient( const bool &isPrevious );

                virtual void setMassChangeVelocityGradientDerivatives( const bool &isPrevious );

                virtual void setMassChangeVelocityGradient( );

                virtual void setPreviousMassChangeVelocityGradient( );

                virtual void setdMassChangeVelocityGradientdDensity( );

                virtual void setdMassChangeVelocityGradientdMassChangeRate( );

                virtual void setdMassChangeVelocityGradientdMassChangeRateGradient( );

                virtual void setdPreviousMassChangeVelocityGradientdPreviousDensity( );

                virtual void setdPreviousMassChangeVelocityGradientdPreviousMassChangeRate( );

                virtual void setdPreviousMassChangeVelocityGradientdPreviousMassChangeRateGradient( );

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

                virtual void setdMassChangeIntermediateVelocityGradientdDensity( );

                virtual void setdMassChangeIntermediateVelocityGradientdMassChangeRate( );

                virtual void setdMassChangeIntermediateVelocityGradientdMassChangeRateGradient( );

                virtual void setdMassChangeIntermediateVelocityGradientdDeformationGradient( );

                virtual void setdMassChangeIntermediateVelocityGradientdSubDeformationGradients( );

                virtual void setPreviousMassChangeIntermediateVelocityGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRateGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( );

                virtual void setdPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( );

                virtual void setMassChangeDeformationGradient( );

                virtual void setMassChangeDeformationGradientDerivatives( const bool &computePrevious );

                virtual void setdMassChangeDeformationGradientdDensity( );

                virtual void setdMassChangeDeformationGradientdMassChangeRate( );

                virtual void setdMassChangeDeformationGradientdMassChangeRateGradient( );

                virtual void setdMassChangeDeformationGradientdDeformationGradient( );

                virtual void setdMassChangeDeformationGradientdSubDeformationGradients( );

                virtual void setdMassChangeDeformationGradientdPreviousDensity( );

                virtual void setdMassChangeDeformationGradientdPreviousMassChangeRate( );

                virtual void setdMassChangeDeformationGradientdPreviousMassChangeRateGradient( );

                virtual void setdMassChangeDeformationGradientdPreviousDeformationGradient( );

                virtual void setdMassChangeDeformationGradientdPreviousSubDeformationGradients( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdAdditionalDOF( ) override;

            private:

                // Friend classes
                friend class tardigradeHydra::massChange::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBase::residualBase;

                using tardigradeHydra::residualBase::setResidual;

                using tardigradeHydra::residualBase::setJacobian;

                using tardigradeHydra::residualBase::setdRdF;

                using tardigradeHydra::residualBase::setdRdT;

                using tardigradeHydra::residualBase::setdRdAdditionalDOF;

                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                unsigned int _massChangeConfigurationIndex;

                floatType _integrationParameter;

                virtual void decomposeParameters( const floatVector &parameters );

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, density,                                         floatType,   unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, previousDensity,                                 floatType,   unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, massChangeRate,                                  floatType,   unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, previousMassChangeRate,                          floatType,   unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, massChangeRateGradient,                          floatVector, unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, previousMassChangeRateGradient,                  floatVector, unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, massDirectionMixingParameter,                    floatType,   unexpectedError                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, massChangeVelocityGradientTrace,                 floatType,   setMassChangeVelocityGradientTrace                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeVelocityGradientTracedDensity,        floatType,   setdMassChangeVelocityGradientTracedDensity        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeVelocityGradientTracedMassChangeRate, floatType,   setdMassChangeVelocityGradientTracedMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMassChangeVelocityGradientTrace,         floatType,   setPreviousMassChangeVelocityGradientTrace         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeVelocityGradientTracedPreviousDensity, floatType, setdPreviousMassChangeVelocityGradientTracedPreviousDensity        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate, floatType, setdPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, unitDirectionVector,                                     floatVector, setUnitDirectionVector                                     ) 

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousUnitDirectionVector,                             floatVector, setPreviousUnitDirectionVector                             ) 

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUnitDirectionVectordMassChangeRateGradient,                 floatVector, setdUnitDirectionVectordMassChangeRateGradient                 ) 

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousUnitDirectionVectordPreviousMassChangeRateGradient, floatVector, setdPreviousUnitDirectionVectordPreviousMassChangeRateGradient ) 

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, massChangeVelocityGradient,                          floatVector, setMassChangeVelocityGradient                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeVelocityGradientdDensity,                 floatVector, setdMassChangeVelocityGradientdDensity                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeVelocityGradientdMassChangeRate,          floatVector, setdMassChangeVelocityGradientdMassChangeRate          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeVelocityGradientdMassChangeRateGradient,  floatVector, setdMassChangeVelocityGradientdMassChangeRateGradient  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMassChangeVelocityGradient,                  floatVector, setPreviousMassChangeVelocityGradient                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeVelocityGradientdPreviousDensity,                floatVector, setdPreviousMassChangeVelocityGradientdPreviousDensity                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeVelocityGradientdPreviousMassChangeRate,         floatVector, setdPreviousMassChangeVelocityGradientdPreviousMassChangeRate          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeVelocityGradientdPreviousMassChangeRateGradient, floatVector, setdPreviousMassChangeVelocityGradientdPreviousMassChangeRateGradient  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingDeformationGradient,                                       floatVector, setPrecedingDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdDeformationGradient,                  floatVector, setdPrecedingDeformationGradientdDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdSubDeformationGradients,              floatVector, setdPrecedingDeformationGradientdSubDeformationGradients )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousPrecedingDeformationGradient,                                       floatVector, setPreviousPrecedingDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, dPreviousPrecedingDeformationGradientdPreviousDeformationGradient,          floatVector, setdPreviousPrecedingDeformationGradientdPreviousDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients,      floatVector, setdPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, massChangeIntermediateVelocityGradient,                             floatVector, setMassChangeIntermediateVelocityGradient         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeIntermediateVelocityGradientdDensity,                    floatVector, setdMassChangeIntermediateVelocityGradientdDensity )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeIntermediateVelocityGradientdMassChangeRate,             floatVector, setdMassChangeIntermediateVelocityGradientdMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeIntermediateVelocityGradientdMassChangeRateGradient,     floatVector, setdMassChangeIntermediateVelocityGradientdMassChangeRateGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeIntermediateVelocityGradientdDeformationGradient,        floatVector, setdMassChangeIntermediateVelocityGradientdDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeIntermediateVelocityGradientdSubDeformationGradients,    floatVector, setdMassChangeIntermediateVelocityGradientdSubDeformationGradients )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMassChangeIntermediateVelocityGradient,                     floatVector, setPreviousMassChangeIntermediateVelocityGradient )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity,                    floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousDensity )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate,             floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRateGradient,     floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRateGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient,        floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients,    floatVector, setdPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, massChangeDeformationGradient,                                                      floatVector, setMassChangeDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeDeformationGradientdDensity,                                             floatVector, setdMassChangeDeformationGradientdDensity )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeDeformationGradientdMassChangeRate,                                      floatVector, setdMassChangeDeformationGradientdMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeDeformationGradientdMassChangeRateGradient,                              floatVector, setdMassChangeDeformationGradientdMassChangeRateGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeDeformationGradientdDeformationGradient,                                 floatVector, setdMassChangeDeformationGradientdDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMassChangeDeformationGradientdSubDeformationGradients,                             floatVector, setdMassChangeDeformationGradientdSubDeformationGradients )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dMassChangeDeformationGradientdPreviousDensity,                                     floatVector, setdMassChangeDeformationGradientdPreviousDensity )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dMassChangeDeformationGradientdPreviousMassChangeRate,                              floatVector, setdMassChangeDeformationGradientdPreviousMassChangeRate )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dMassChangeDeformationGradientdPreviousMassChangeRateGradient,                      floatVector, setdMassChangeDeformationGradientdPreviousMassChangeRateGradient )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dMassChangeDeformationGradientdPreviousDeformationGradient,                         floatVector, setdMassChangeDeformationGradientdPreviousDeformationGradient )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dMassChangeDeformationGradientdPreviousSubDeformationGradients,                     floatVector, setdMassChangeDeformationGradientdPreviousSubDeformationGradients )

        };

    }

}

#endif
