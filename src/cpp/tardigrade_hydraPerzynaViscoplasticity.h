/**
  ******************************************************************************
  * \file tardigrade_hydraPerzynaViscoplasticity.h
  ******************************************************************************
  * An implementation of perzynaViscoplasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_PERZYNA_VISCOPLASTICITY_H
#define TARDIGRADE_HYDRA_PERZYNA_VISCOPLASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace perzynaViscoplasticity{

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

                    _parameters = parameters;

                }

            public:

                 // Friend classes
                friend class tardigradeHydra::perzynaViscoplasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

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

                virtual void setdYieldFunctiondStateVariables( );

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

                virtual void setdPreviousYieldFunctiondPreviousStateVariables( );

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

                virtual void setdYieldFunctiondStateVariables( const bool isPrevious );

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

                //! Get a reference to the parameter vector
                const floatVector * getParameters( ){ return &_parameters; }

                //! Set the Perzyna parameters
                virtual void setPerzynaParameters( ){ decomposeParameters( *getParameters( ) ); }

                //! Set the drag stress parameters
                virtual void setDragStressParameters( ){ decomposeParameters( *getParameters( ) ); }

                //! Set the thermal parameters
                virtual void setThermalParameters( ){ decomposeParameters( *getParameters( ) ); }

                //! Set the yield parameters
                virtual void setYieldParameters( ){ decomposeParameters( *getParameters( ) ); }

                //! Set the flow parameters
                virtual void setFlowParameters( ){ decomposeParameters( *getParameters( ) ); }

                //! Set the hardening parameters
                virtual void setHardeningParameters( ){ decomposeParameters( *getParameters( ) ); }

            private:

                unsigned int _plasticConfigurationIndex;

                std::vector< unsigned int > _stateVariableIndices;

                floatType _integrationParameter;

                floatVector _parameters;

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, drivingStress,                                               secondOrderTensor, setDrivingStress                                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDrivingStress,                                       secondOrderTensor, setPreviousDrivingStress                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdCauchyStress,                                 fourthOrderTensor, setdDrivingStressdCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdF,                                            fourthOrderTensor, setdDrivingStressdF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdSubFs,                                        floatVector, setdDrivingStressdSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousCauchyStress,                 fourthOrderTensor, setdPreviousDrivingStressdPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousF,                            fourthOrderTensor, setdPreviousDrivingStressdPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousSubFs,                        floatVector, setdPreviousDrivingStressdPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, flowDirection,                                               secondOrderTensor, setFlowDirection                                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousFlowDirection,                                       secondOrderTensor, setPreviousFlowDirection                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondCauchyStress,                                 fourthOrderTensor, setdFlowDirectiondCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondF,                                            fourthOrderTensor, setdFlowDirectiondF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondSubFs,                                        floatVector, setdFlowDirectiondSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousCauchyStress,                 fourthOrderTensor, setdPreviousFlowDirectiondPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousF,                            fourthOrderTensor, setdPreviousFlowDirectiondPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousSubFs,                        floatVector, setdPreviousFlowDirectiondPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, yieldFunction,                                               floatType,   setYieldFunction                                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondCauchyStress,                                 secondOrderTensor, setdYieldFunctiondCauchyStress                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondF,                                            secondOrderTensor, setdYieldFunctiondF                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondSubFs,                                        floatVector, setdYieldFunctiondSubFs                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondStateVariables,                               floatVector, setdYieldFunctiondStateVariables                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousYieldFunction,                                       floatType,   setPreviousYieldFunction                                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousCauchyStress,                 secondOrderTensor, setdPreviousYieldFunctiondPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousF,                            secondOrderTensor, setdPreviousYieldFunctiondPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousSubFs,                        floatVector, setdPreviousYieldFunctiondPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousStateVariables,               floatVector, setdPreviousYieldFunctiondPreviousStateVariables               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticThermalMultiplier,                                    floatType,   setPlasticThermalMultiplier                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticThermalMultiplierdT,                                 floatType,   setdPlasticThermalMultiplierdT                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticThermalMultiplier,                            floatType,   setPreviousPlasticThermalMultiplier                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticThermalMultiplierdPreviousT,                 floatType,   setdPreviousPlasticThermalMultiplierdPreviousT                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dragStress,                                                  floatType,   setDragStress                                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDragStressdStateVariables,                                  floatVector, setdDragStressdStateVariables                                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDragStress,                                          floatType,   setPreviousDragStress                                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDragStressdPreviousStateVariables,                  floatVector, setdPreviousDragStressdPreviousStateVariables                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, hardeningFunction,                                           floatVector, setHardeningFunction                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHardeningFunctiondStateVariables,                           floatVector, setdHardeningFunctiondStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHardeningFunction,                                   floatVector, setPreviousHardeningFunction                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousHardeningFunctiondPreviousStateVariables,           floatVector, setdPreviousHardeningFunctiondPreviousStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMultiplier,                                           floatType,   setPlasticMultiplier                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdCauchyStress,                             secondOrderTensor, setdPlasticMultiplierdCauchyStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdF,                                        secondOrderTensor, setdPlasticMultiplierdF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdSubFs,                                    floatVector, setdPlasticMultiplierdSubFs                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdT,                                        floatType,   setdPlasticMultiplierdT                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdStateVariables,                           floatVector, setdPlasticMultiplierdStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMultiplier,                                   floatType,   setPreviousPlasticMultiplier                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousCauchyStress,             secondOrderTensor, setdPreviousPlasticMultiplierdPreviousCauchyStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousF,                        secondOrderTensor, setdPreviousPlasticMultiplierdPreviousF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousSubFs,                    floatVector, setdPreviousPlasticMultiplierdPreviousSubFs                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousT,                        floatType,   setdPreviousPlasticMultiplierdPreviousT                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousStateVariables,           floatVector, setdPreviousPlasticMultiplierdPreviousStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, velocityGradient,                                            secondOrderTensor, setVelocityGradient                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdCauchyStress,                              fourthOrderTensor, setdVelocityGradientdCauchyStress                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdF,                                         fourthOrderTensor, setdVelocityGradientdF                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdSubFs,                                     floatVector, setdVelocityGradientdSubFs                                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdT,                                         secondOrderTensor, setdVelocityGradientdT                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdStateVariables,                            floatVector, setdVelocityGradientdStateVariables                            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousVelocityGradient,                                    secondOrderTensor, setPreviousVelocityGradient                                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousCauchyStress,              fourthOrderTensor, setdPreviousVelocityGradientdPreviousCauchyStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousF,                         fourthOrderTensor, setdPreviousVelocityGradientdPreviousF                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousSubFs,                     floatVector, setdPreviousVelocityGradientdPreviousSubFs                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousT,                         secondOrderTensor, setdPreviousVelocityGradientdPreviousT                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousStateVariables,            floatVector, setdPreviousVelocityGradientdPreviousStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableEvolutionRates,                                 floatVector, setStateVariableEvolutionRates                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdCauchyStress,                   floatVector, setdStateVariableEvolutionRatesdCauchyStress                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdF,                              floatVector, setdStateVariableEvolutionRatesdF                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdSubFs,                          floatVector, setdStateVariableEvolutionRatesdSubFs                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdT,                              floatVector, setdStateVariableEvolutionRatesdT                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdStateVariables,                 floatVector, setdStateVariableEvolutionRatesdStateVariables                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariableEvolutionRates,                         floatVector, setPreviousStateVariableEvolutionRates                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousCauchyStress,   floatVector, setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousF,              floatVector, setdPreviousStateVariableEvolutionRatesdPreviousF              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousSubFs,          floatVector, setdPreviousStateVariableEvolutionRatesdPreviousSubFs          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousT,              floatVector, setdPreviousStateVariableEvolutionRatesdPreviousT              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousStateVariables, floatVector, setdPreviousStateVariableEvolutionRatesdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticDeformationGradient,                                  secondOrderTensor, setPlasticDeformationGradient                                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdCauchyStress,                    fourthOrderTensor, setdPlasticDeformationGradientdCauchyStress                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdF,                               fourthOrderTensor, setdPlasticDeformationGradientdF                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdSubFs,                           floatVector, setdPlasticDeformationGradientdSubFs                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdT,                               secondOrderTensor, setdPlasticDeformationGradientdT                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdStateVariables,                  floatVector, setdPlasticDeformationGradientdStateVariables                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousCauchyStress,            fourthOrderTensor, setdPlasticDeformationGradientdPreviousCauchyStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousF,                       fourthOrderTensor, setdPlasticDeformationGradientdPreviousF                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousSubFs,                   floatVector, setdPlasticDeformationGradientdPreviousSubFs                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousT,                       secondOrderTensor, setdPlasticDeformationGradientdPreviousT                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousStateVariables,          floatVector, setdPlasticDeformationGradientdPreviousStateVariables          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStateVariables,                                       floatVector, setPlasticStateVariables                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdCauchyStress,                         floatVector, setdPlasticStateVariablesdCauchyStress                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdF,                                    floatVector, setdPlasticStateVariablesdF                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdSubFs,                                floatVector, setdPlasticStateVariablesdSubFs                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdT,                                    floatVector, setdPlasticStateVariablesdT                                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdStateVariables,                       floatVector, setdPlasticStateVariablesdStateVariables                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousCauchyStress,                 floatVector, setdPlasticStateVariablesdPreviousCauchyStress                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousF,                            floatVector, setdPlasticStateVariablesdPreviousF                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousSubFs,                        floatVector, setdPlasticStateVariablesdPreviousSubFs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousT,                            floatVector, setdPlasticStateVariablesdPreviousT                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousStateVariables,               floatVector, setdPlasticStateVariablesdPreviousStateVariables               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariables,                                              floatVector, setStateVariables                                              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariables,                                      floatVector, setPreviousStateVariables                                      )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, perzynaParameters,                                           floatVector, setPerzynaParameters                                           )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, dragStressParameters,                                        floatVector, setDragStressParameters                                        )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, thermalParameters,                                           floatVector, setThermalParameters                                           )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, yieldParameters,                                             floatVector, setYieldParameters                                             )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, flowParameters,                                              floatVector, setFlowParameters                                              )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, hardeningParameters,                                         floatVector, setHardeningParameters                                         )

        };

    }

}

#endif
