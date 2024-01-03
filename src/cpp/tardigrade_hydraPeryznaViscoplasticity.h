/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaViscoplasticity.h
  ******************************************************************************
  * An implementation of peryznaViscoplasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_PERYZNA_VISCOPLASTICITY_H
#define TARDIGRADE_HYDRA_PERYZNA_VISCOPLASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace peryznaViscoplasticity{

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
         * A class which defines a viscoplastic residual
         */
        class residual : public tardigradeHydra::residualBase{

            public:

                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const floatType integrationParameter = 0.5 ) : tardigradeHydra::residualBase( hydra, numEquations ){
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
                     * \param &integrationParameter: The integration parameter for the function. 0 is explicit, 1 is implicit.
                     */

                    _plasticConfigurationIndex = plasticConfigurationIndex;

                    _stateVariableIndices = stateVariableIndices;

                    _integrationParameter = integrationParameter;

                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameters( parameters ) );

                }

            public:

                 // Friend classes
                friend class tardigradeHydra::peryznaViscoplasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

                using tardigradeHydra::residualBase::residualBase;

                using tardigradeHydra::residualBase::setResidual;

                using tardigradeHydra::residualBase::setJacobian;

                using tardigradeHydra::residualBase::setdRdF;

                using tardigradeHydra::residualBase::setdRdT;

                using tardigradeHydra::residualBase::setAdditionalDerivatives;

                const unsigned int* getPlasticConfigurationIndex( );

                const std::vector< unsigned int >* getStateVariableIndices( );

                const floatType* getIntegrationParameter( );

            protected:

                virtual void setDrivingStress( );

                virtual void setdDrivingStressdCauchyStress( );

                virtual void setdDrivingStressdF( );

                virtual void setdDrivingStressdSubFs( );

                virtual void setFlowDirection( );

                virtual void setdFlowDirectiondCauchyStress( );

                virtual void setdFlowDirectiondF( );

                virtual void setdFlowDirectiondSubFs( );

                virtual void setYieldFunction( );

                virtual void setdYieldFunctiondCauchyStress( );

                virtual void setdYieldFunctiondF( );

                virtual void setdYieldFunctiondSubFs( );

                virtual void setPlasticThermalMultiplier( );

                virtual void setdPlasticThermalMultiplierdT( );

                virtual void setDragStress( );

                virtual void setdDragStressdStateVariables( );

                virtual void setHardeningFunction( );

                virtual void setdHardeningFunctiondStateVariables( );

                virtual void setPlasticMultiplier( );

                virtual void setdPlasticMultiplierdCauchyStress( );

                virtual void setdPlasticMultiplierdF( );

                virtual void setdPlasticMultiplierdSubFs( );

                virtual void setdPlasticMultiplierdT( );

                virtual void setdPlasticMultiplierdStateVariables( );

                virtual void setVelocityGradient( );

                virtual void setdVelocityGradientdCauchyStress( );

                virtual void setdVelocityGradientdF( );

                virtual void setdVelocityGradientdSubFs( );

                virtual void setdVelocityGradientdT( );

                virtual void setdVelocityGradientdStateVariables( );

                virtual void setStateVariableEvolutionRates( );

                virtual void setdStateVariableEvolutionRatesdCauchyStress( );

                virtual void setdStateVariableEvolutionRatesdF( );

                virtual void setdStateVariableEvolutionRatesdSubFs( );

                virtual void setdStateVariableEvolutionRatesdT( );

                virtual void setdStateVariableEvolutionRatesdStateVariables( );

                virtual void setPreviousDrivingStress( );

                virtual void setdPreviousDrivingStressdPreviousCauchyStress( );

                virtual void setdPreviousDrivingStressdPreviousF( );

                virtual void setdPreviousDrivingStressdPreviousSubFs( );

                virtual void setPreviousFlowDirection( );

                virtual void setdPreviousFlowDirectiondPreviousCauchyStress( );

                virtual void setdPreviousFlowDirectiondPreviousF( );

                virtual void setdPreviousFlowDirectiondPreviousSubFs( );

                virtual void setPreviousYieldFunction( );

                virtual void setdPreviousYieldFunctiondPreviousCauchyStress( );

                virtual void setdPreviousYieldFunctiondPreviousF( );

                virtual void setdPreviousYieldFunctiondPreviousSubFs( );

                virtual void setPreviousPlasticThermalMultiplier( );

                virtual void setdPreviousPlasticThermalMultiplierdPreviousT( );

                virtual void setPreviousDragStress( );

                virtual void setdPreviousDragStressdPreviousStateVariables( );

                virtual void setPreviousHardeningFunction( );

                virtual void setdPreviousHardeningFunctiondPreviousStateVariables( );

                virtual void setPreviousPlasticMultiplier( );

                virtual void setdPreviousPlasticMultiplierdPreviousCauchyStress( );

                virtual void setdPreviousPlasticMultiplierdPreviousF( );

                virtual void setdPreviousPlasticMultiplierdPreviousSubFs( );

                virtual void setdPreviousPlasticMultiplierdPreviousT( );

                virtual void setdPreviousPlasticMultiplierdPreviousStateVariables( );

                virtual void setPreviousVelocityGradient( );

                virtual void setdPreviousVelocityGradientdPreviousCauchyStress( );

                virtual void setdPreviousVelocityGradientdPreviousF( );

                virtual void setdPreviousVelocityGradientdPreviousSubFs( );

                virtual void setdPreviousVelocityGradientdPreviousT( );

                virtual void setdPreviousVelocityGradientdPreviousStateVariables( );

                virtual void setPreviousStateVariableEvolutionRates( );

                virtual void setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

                virtual void setdPreviousStateVariableEvolutionRatesdPreviousF( );

                virtual void setdPreviousStateVariableEvolutionRatesdPreviousSubFs( );

                virtual void setdPreviousStateVariableEvolutionRatesdPreviousT( );

                virtual void setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( );

                virtual void setDrivingStress( const bool isPrevious );

                virtual void setDrivingStressDerivatives( const bool isPrevious );

                virtual void setdDrivingStressdCauchyStress( const bool isPrevious );

                virtual void setdDrivingStressdF( const bool isPrevious );

                virtual void setdDrivingStressdSubFs( const bool isPrevious );

                virtual void setFlowDirection( const bool isPrevious );

                virtual void setFlowDirectionDerivatives( const bool isPrevious );

                virtual void setdFlowDirectiondCauchyStress( const bool isPrevious );

                virtual void setdFlowDirectiondF( const bool isPrevious );

                virtual void setdFlowDirectiondSubFs( const bool isPrevious );

                virtual void setYieldFunction( const bool isPrevious);

                virtual void setYieldFunctionDerivatives( const bool isPrevious );

                virtual void setdYieldFunctiondCauchyStress( const bool isPrevious );

                virtual void setdYieldFunctiondF( const bool isPrevious );

                virtual void setdYieldFunctiondSubFs( const bool isPrevious );

                virtual void setPlasticThermalMultiplier( const bool isPrevious );

                virtual void setPlasticThermalMultiplierDerivatives( const bool isPrevious );

                virtual void setdPlasticThermalMultiplierdT( const bool isPrevious );

                virtual void setDragStress( const bool isPrevious );

                virtual void setDragStressDerivatives( const bool isPrevious );

                virtual void setdDragStressdStateVariables( const bool isPrevious );

                virtual void setHardeningFunction( const bool isPrevious );

                virtual void setHardeningFunctionDerivatives( const bool isPrevious );

                virtual void setdHardeningFunctiondStateVariables( const bool isPrevious );

                virtual void setPlasticMultiplier( const bool isPrevious );

                virtual void setPlasticMultiplierDerivatives( const bool isPrevious );

                virtual void setdPlasticMultiplierdCauchyStress( const bool isPrevious );

                virtual void setdPlasticMultiplierdF( const bool isPrevious );

                virtual void setdPlasticMultiplierdSubFs( const bool isPrevious );

                virtual void setdPlasticMultiplierdT( const bool isPrevious );

                virtual void setdPlasticMultiplierdStateVariables( const bool isPrevious );

                virtual void setVelocityGradient( const bool isPrevious );

                virtual void setVelocityGradientDerivatives( const bool isPrevious );

                virtual void setdVelocityGradientdCauchyStress( const bool isPrevious );

                virtual void setdVelocityGradientdF( const bool isPrevious );

                virtual void setdVelocityGradientdSubFs( const bool isPrevious );

                virtual void setdVelocityGradientdT( const bool isPrevious );

                virtual void setdVelocityGradientdStateVariables( const bool isPrevious );

                virtual void setStateVariableEvolutionRates( const bool isPrevious );

                virtual void setStateVariableEvolutionRateDerivatives( const bool isPrevious );

                virtual void setdStateVariableEvolutionRatesdCauchyStress( const bool isPrevious );

                virtual void setdStateVariableEvolutionRatesdF( const bool isPrevious );

                virtual void setdStateVariableEvolutionRatesdSubFs( const bool isPrevious );

                virtual void setdStateVariableEvolutionRatesdT( const bool isPrevious );

                virtual void setdStateVariableEvolutionRatesdStateVariables( const bool isPrevious );

                virtual void setPlasticDeformationGradient( );

                virtual void setPlasticDeformationGradientDerivatives( const bool setPreviousDerivatives );

                virtual void setdPlasticDeformationGradientdCauchyStress( );

                virtual void setdPlasticDeformationGradientdF( );

                virtual void setdPlasticDeformationGradientdSubFs( );

                virtual void setdPlasticDeformationGradientdT( );

                virtual void setdPlasticDeformationGradientdStateVariables( );

                virtual void setdPlasticDeformationGradientdPreviousCauchyStress( );

                virtual void setdPlasticDeformationGradientdPreviousF( );

                virtual void setdPlasticDeformationGradientdPreviousSubFs( );

                virtual void setdPlasticDeformationGradientdPreviousT( );

                virtual void setdPlasticDeformationGradientdPreviousStateVariables( );

                virtual void setPlasticStateVariables( );

                virtual void setPlasticStateVariableDerivatives( const bool setPreviousDerivatives );

                virtual void setdPlasticStateVariablesdCauchyStress( );

                virtual void setdPlasticStateVariablesdF( );

                virtual void setdPlasticStateVariablesdSubFs( );

                virtual void setdPlasticStateVariablesdT( );

                virtual void setdPlasticStateVariablesdStateVariables( );

                virtual void setdPlasticStateVariablesdPreviousCauchyStress( );

                virtual void setdPlasticStateVariablesdPreviousF( );

                virtual void setdPlasticStateVariablesdPreviousSubFs( );

                virtual void setdPlasticStateVariablesdPreviousT( );

                virtual void setdPlasticStateVariablesdPreviousStateVariables( );

                virtual void setStateVariables( );

                virtual void setStateVariables( const bool isPrevious );

                virtual void setPreviousStateVariables( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdT( ) override;

                virtual void decomposeParameters( const floatVector &parameters );

            private:

                unsigned int _plasticConfigurationIndex;

                std::vector< unsigned int > _stateVariableIndices;

                floatType _integrationParameter;

                //! Throw an error if the Peryzna parameters haven't been defined
                void peryznaError( ){ throw std::runtime_error( "Peryzna parameters not defined but required" ); }

                //! Throw an error if the drag stress parameters haven't been defined
                void dragStressError( ){ throw std::runtime_error( "Drag stress parameters not defined but required" ); }

                //! Throw an error if the thermal parameters haven't been defined
                void thermalError( ){ throw std::runtime_error( "Thermal parameters not defined but required" ); }

                //! Throw an error if the yield parameters haven't been defined
                void yieldError( ){ throw std::runtime_error( "Yield parameters not defined but required" ); }

                //! Throw an error if the flow parameters haven't been defined
                void flowError( ){ throw std::runtime_error( "Flow parameters not defined but required" ); }

                //! Throw an error if the hardening parameters haven't been defined
                void hardeningError( ){ throw std::runtime_error( "Hardening parameters not defined but required" ); }

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, drivingStress,                                               floatVector, setDrivingStress                                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDrivingStress,                                       floatVector, setPreviousDrivingStress                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdCauchyStress,                                 floatMatrix, setdDrivingStressdCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdF,                                            floatMatrix, setdDrivingStressdF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdSubFs,                                        floatMatrix, setdDrivingStressdSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousCauchyStress,                 floatMatrix, setdPreviousDrivingStressdPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousF,                            floatMatrix, setdPreviousDrivingStressdPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousSubFs,                        floatMatrix, setdPreviousDrivingStressdPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, flowDirection,                                               floatVector, setFlowDirection                                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousFlowDirection,                                       floatVector, setPreviousFlowDirection                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondCauchyStress,                                 floatMatrix, setdFlowDirectiondCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondF,                                            floatMatrix, setdFlowDirectiondF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondSubFs,                                        floatMatrix, setdFlowDirectiondSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousCauchyStress,                 floatMatrix, setdPreviousFlowDirectiondPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousF,                            floatMatrix, setdPreviousFlowDirectiondPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousSubFs,                        floatMatrix, setdPreviousFlowDirectiondPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, yieldFunction,                                               floatType,   setYieldFunction                                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondCauchyStress,                                 floatVector, setdYieldFunctiondCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondF,                                            floatVector, setdYieldFunctiondF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondSubFs,                                        floatVector, setdYieldFunctiondSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousYieldFunction,                                       floatType,   setPreviousYieldFunction                                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousCauchyStress,                 floatVector, setdPreviousYieldFunctiondPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousF,                            floatVector, setdPreviousYieldFunctiondPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousSubFs,                        floatVector, setdPreviousYieldFunctiondPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticThermalMultiplier,                                    floatType,   setPlasticThermalMultiplier                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticThermalMultiplierdT,                                 floatType,   setdPlasticThermalMultiplierdT                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticThermalMultiplier,                            floatType,   setPreviousPlasticThermalMultiplier                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticThermalMultiplierdPreviousT,                 floatType,   setdPreviousPlasticThermalMultiplierdPreviousT                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dragStress,                                                  floatType,   setDragStress                                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDragStressdStateVariables,                                  floatVector, setdDragStressdStateVariables                                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDragStress,                                          floatType,   setPreviousDragStress                                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDragStressdPreviousStateVariables,                  floatVector, setdPreviousDragStressdPreviousStateVariables                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, hardeningFunction,                                           floatType,   setHardeningFunction                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHardeningFunctiondStateVariables,                           floatVector, setdHardeningFunctiondStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHardeningFunction,                                   floatType,   setPreviousHardeningFunction                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousHardeningFunctiondPreviousStateVariables,           floatVector, setdPreviousHardeningFunctiondPreviousStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMultiplier,                                           floatType,   setPlasticMultiplier                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdCauchyStress,                             floatVector, setdPlasticMultiplierdCauchyStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdF,                                        floatVector, setdPlasticMultiplierdF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdSubFs,                                    floatVector, setdPlasticMultiplierdSubFs                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdT,                                        floatType,   setdPlasticMultiplierdT                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdStateVariables,                           floatVector, setdPlasticMultiplierdStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMultiplier,                                   floatType,   setPreviousPlasticMultiplier                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousCauchyStress,             floatVector, setdPreviousPlasticMultiplierdPreviousCauchyStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousF,                        floatVector, setdPreviousPlasticMultiplierdPreviousF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousSubFs,                    floatVector, setdPreviousPlasticMultiplierdPreviousSubFs                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousT,                        floatType,   setdPreviousPlasticMultiplierdPreviousT                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousStateVariables,           floatVector, setdPreviousPlasticMultiplierdPreviousStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, velocityGradient,                                            floatVector, setVelocityGradient                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdCauchyStress,                              floatMatrix, setdVelocityGradientdCauchyStress                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdF,                                         floatMatrix, setdVelocityGradientdF                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdSubFs,                                     floatMatrix, setdVelocityGradientdSubFs                                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdT,                                         floatVector, setdVelocityGradientdT                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdStateVariables,                            floatMatrix, setdVelocityGradientdStateVariables                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousVelocityGradient,                                    floatVector, setPreviousVelocityGradient                                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousCauchyStress,              floatMatrix, setdPreviousVelocityGradientdPreviousCauchyStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousF,                         floatMatrix, setdPreviousVelocityGradientdPreviousF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousSubFs,                     floatMatrix, setdPreviousVelocityGradientdPreviousSubFs                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousT,                         floatVector, setdPreviousVelocityGradientdPreviousT                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousStateVariables,            floatMatrix, setdPreviousVelocityGradientdPreviousStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableEvolutionRates,                                 floatVector, setStateVariableEvolutionRates                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdCauchyStress,                   floatMatrix, setdStateVariableEvolutionRatesdCauchyStress                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdF,                              floatMatrix, setdStateVariableEvolutionRatesdF                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdSubFs,                          floatMatrix, setdStateVariableEvolutionRatesdSubFs                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdT,                              floatVector, setdStateVariableEvolutionRatesdT                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdStateVariables,                 floatMatrix, setdStateVariableEvolutionRatesdStateVariables                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariableEvolutionRates,                         floatVector, setPreviousStateVariableEvolutionRates                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousCauchyStress,   floatMatrix, setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousF,              floatMatrix, setdPreviousStateVariableEvolutionRatesdPreviousF              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousSubFs,          floatMatrix, setdPreviousStateVariableEvolutionRatesdPreviousSubFs          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousT,              floatVector, setdPreviousStateVariableEvolutionRatesdPreviousT              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousStateVariables, floatMatrix, setdPreviousStateVariableEvolutionRatesdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticDeformationGradient,                                  floatVector, setPlasticDeformationGradient                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdCauchyStress,                    floatMatrix, setdPlasticDeformationGradientdCauchyStress                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdF,                               floatMatrix, setdPlasticDeformationGradientdF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdSubFs,                           floatMatrix, setdPlasticDeformationGradientdSubFs                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdT,                               floatVector, setdPlasticDeformationGradientdT                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdStateVariables,                  floatMatrix, setdPlasticDeformationGradientdStateVariables                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousCauchyStress,            floatMatrix, setdPlasticDeformationGradientdPreviousCauchyStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousF,                       floatMatrix, setdPlasticDeformationGradientdPreviousF                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousSubFs,                   floatMatrix, setdPlasticDeformationGradientdPreviousSubFs                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousT,                       floatVector, setdPlasticDeformationGradientdPreviousT                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousStateVariables,          floatMatrix, setdPlasticDeformationGradientdPreviousStateVariables          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStateVariables,                                       floatVector, setPlasticStateVariables                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdCauchyStress,                         floatMatrix, setdPlasticStateVariablesdCauchyStress                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdF,                                    floatMatrix, setdPlasticStateVariablesdF                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdSubFs,                                floatMatrix, setdPlasticStateVariablesdSubFs                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdT,                                    floatVector, setdPlasticStateVariablesdT                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdStateVariables,                       floatMatrix, setdPlasticStateVariablesdStateVariables                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousCauchyStress,                 floatMatrix, setdPlasticStateVariablesdPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousF,                            floatMatrix, setdPlasticStateVariablesdPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousSubFs,                        floatMatrix, setdPlasticStateVariablesdPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousT,                            floatVector, setdPlasticStateVariablesdPreviousT                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousStateVariables,               floatMatrix, setdPlasticStateVariablesdPreviousStateVariables               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariables,                                              floatVector, setStateVariables                                              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariables,                                      floatVector, setPreviousStateVariables                                      )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, peryznaParameters,                                           floatVector, peryznaError                                                   )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, dragStressParameters,                                        floatVector, dragStressError                                                )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, thermalParameters,                                           floatVector, thermalError                                                   )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, yieldParameters,                                             floatVector, yieldError                                                     )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, flowParameters,                                              floatVector, flowError                                                      )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, hardeningParameters,                                         floatVector, hardeningError                                                 )

        };

    }

}

#endif
