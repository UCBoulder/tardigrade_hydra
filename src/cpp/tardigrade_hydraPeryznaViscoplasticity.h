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

                const floatVector* getDrivingStress( );

                const floatMatrix* getdDrivingStressdCauchyStress( );

                const floatMatrix* getdDrivingStressdF( );

                const floatMatrix* getdDrivingStressdSubFs( );

                const floatVector* getFlowDirection( );

                const floatMatrix* getdFlowDirectiondCauchyStress( );

                const floatMatrix* getdFlowDirectiondF( );

                const floatMatrix* getdFlowDirectiondSubFs( );

                const floatType* getYieldFunction( );

                const floatVector* getdYieldFunctiondCauchyStress( );

                const floatVector* getdYieldFunctiondF( );

                const floatVector* getdYieldFunctiondSubFs( );

                const floatType* getPlasticThermalMultiplier( );

                const floatType* getdPlasticThermalMultiplierdT( );

                const floatType* getDragStress( );

                const floatVector* getdDragStressdStateVariables( );

                const floatType* getHardeningFunction( );

                const floatVector* getdHardeningFunctiondStateVariables( );

                const floatType* getPlasticMultiplier( );

                const floatVector* getdPlasticMultiplierdCauchyStress( );

                const floatVector* getdPlasticMultiplierdF( );

                const floatVector* getdPlasticMultiplierdSubFs( );

                const floatType* getdPlasticMultiplierdT( );

                const floatVector* getdPlasticMultiplierdStateVariables( );

                const floatVector* getVelocityGradient( );

                const floatMatrix* getdVelocityGradientdCauchyStress( );

                const floatMatrix* getdVelocityGradientdF( );

                const floatMatrix* getdVelocityGradientdSubFs( );

                const floatVector* getdVelocityGradientdT( );

                const floatMatrix* getdVelocityGradientdStateVariables( );

                const floatVector* getStateVariableEvolutionRates( );

                const floatMatrix* getdStateVariableEvolutionRatesdCauchyStress( );

                const floatMatrix* getdStateVariableEvolutionRatesdF( );

                const floatMatrix* getdStateVariableEvolutionRatesdSubFs( );

                const floatVector* getdStateVariableEvolutionRatesdT( );

                const floatMatrix* getdStateVariableEvolutionRatesdStateVariables( );

                const floatVector* getPreviousDrivingStress( );

                const floatMatrix* getdPreviousDrivingStressdPreviousCauchyStress( );

                const floatMatrix* getdPreviousDrivingStressdPreviousF( );

                const floatMatrix* getdPreviousDrivingStressdPreviousSubFs( );

                const floatVector* getPreviousFlowDirection( );

                const floatMatrix* getdPreviousFlowDirectiondPreviousCauchyStress( );

                const floatMatrix* getdPreviousFlowDirectiondPreviousF( );

                const floatMatrix* getdPreviousFlowDirectiondPreviousSubFs( );

                const floatType* getPreviousYieldFunction( );

                const floatVector* getdPreviousYieldFunctiondPreviousCauchyStress( );

                const floatVector* getdPreviousYieldFunctiondPreviousF( );

                const floatVector* getdPreviousYieldFunctiondPreviousSubFs( );

                const floatType* getPreviousPlasticThermalMultiplier( );

                const floatType* getdPreviousPlasticThermalMultiplierdPreviousT( );

                const floatType* getPreviousDragStress( );

                const floatVector* getdPreviousDragStressdPreviousStateVariables( );

                const floatType* getPreviousHardeningFunction( );

                const floatVector* getdPreviousHardeningFunctiondPreviousStateVariables( );

                const floatType* getPreviousPlasticMultiplier( );

                const floatVector* getdPreviousPlasticMultiplierdPreviousCauchyStress( );

                const floatVector* getdPreviousPlasticMultiplierdPreviousF( );

                const floatVector* getdPreviousPlasticMultiplierdPreviousSubFs( );

                const floatType* getdPreviousPlasticMultiplierdPreviousT( );

                const floatVector* getdPreviousPlasticMultiplierdPreviousStateVariables( );

                const floatVector* getPreviousVelocityGradient( );

                const floatMatrix* getdPreviousVelocityGradientdPreviousCauchyStress( );

                const floatMatrix* getdPreviousVelocityGradientdPreviousF( );

                const floatMatrix* getdPreviousVelocityGradientdPreviousSubFs( );

                const floatVector* getdPreviousVelocityGradientdPreviousT( );

                const floatMatrix* getdPreviousVelocityGradientdPreviousStateVariables( );

                const floatVector* getPreviousStateVariableEvolutionRates( );

                const floatMatrix* getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

                const floatMatrix* getdPreviousStateVariableEvolutionRatesdPreviousF( );

                const floatMatrix* getdPreviousStateVariableEvolutionRatesdPreviousSubFs( );

                const floatVector* getdPreviousStateVariableEvolutionRatesdPreviousT( );

                const floatMatrix* getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( );

                const floatVector* getPlasticDeformationGradient( );

                const floatMatrix* getdPlasticDeformationGradientdCauchyStress( );

                const floatMatrix* getdPlasticDeformationGradientdF( );

                const floatMatrix* getdPlasticDeformationGradientdSubFs( );

                const floatVector* getdPlasticDeformationGradientdT( );

                const floatMatrix* getdPlasticDeformationGradientdStateVariables( );

                const floatMatrix* getdPlasticDeformationGradientdPreviousCauchyStress( );

                const floatMatrix* getdPlasticDeformationGradientdPreviousF( );

                const floatMatrix* getdPlasticDeformationGradientdPreviousSubFs( );

                const floatVector* getdPlasticDeformationGradientdPreviousT( );

                const floatMatrix* getdPlasticDeformationGradientdPreviousStateVariables( );

                const floatVector* getPlasticStateVariables( );

                const floatMatrix* getdPlasticStateVariablesdCauchyStress( );

                const floatMatrix* getdPlasticStateVariablesdF( );

                const floatMatrix* getdPlasticStateVariablesdSubFs( );

                const floatVector* getdPlasticStateVariablesdT( );

                const floatMatrix* getdPlasticStateVariablesdStateVariables( );

                const floatMatrix* getdPlasticStateVariablesdPreviousCauchyStress( );

                const floatMatrix* getdPlasticStateVariablesdPreviousF( );

                const floatMatrix* getdPlasticStateVariablesdPreviousSubFs( );

                const floatVector* getdPlasticStateVariablesdPreviousT( );

                const floatMatrix* getdPlasticStateVariablesdPreviousStateVariables( );

                const floatVector* getStateVariables( );

                const floatVector* getPreviousStateVariables( );

                const floatVector* getPeryznaParameters( );

                const floatVector* getDragStressParameters( );

                const floatVector* getThermalParameters( );

                const floatVector* getYieldParameters( );

                const floatVector* getFlowParameters( );

                const floatVector* getHardeningParameters( );

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

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, drivingStress,                                               floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDrivingStress,                                       floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdCauchyStress,                                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdF,                                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDrivingStressdSubFs,                                        floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousCauchyStress,                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousF,                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDrivingStressdPreviousSubFs,                        floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, flowDirection,                                               floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousFlowDirection,                                       floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondCauchyStress,                                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondF,                                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFlowDirectiondSubFs,                                        floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousCauchyStress,                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousF,                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousFlowDirectiondPreviousSubFs,                        floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, yieldFunction,                                               floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondCauchyStress,                                 floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondF,                                            floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dYieldFunctiondSubFs,                                        floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousYieldFunction,                                       floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousCauchyStress,                 floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousF,                            floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousYieldFunctiondPreviousSubFs,                        floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticThermalMultiplier,                                    floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticThermalMultiplierdT,                                 floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticThermalMultiplier,                            floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticThermalMultiplierdPreviousT,                 floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dragStress,                                                  floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dDragStressdStateVariables,                                  floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousDragStress,                                          floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousDragStressdPreviousStateVariables,                  floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, hardeningFunction,                                           floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHardeningFunctiondStateVariables,                           floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHardeningFunction,                                   floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousHardeningFunctiondPreviousStateVariables,           floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMultiplier,                                           floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdCauchyStress,                             floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdF,                                        floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdSubFs,                                    floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdT,                                        floatType   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMultiplierdStateVariables,                           floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMultiplier,                                   floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousCauchyStress,             floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousF,                        floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousSubFs,                    floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousT,                        floatType   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousPlasticMultiplierdPreviousStateVariables,           floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, velocityGradient,                                            floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdCauchyStress,                              floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdF,                                         floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdSubFs,                                     floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdT,                                         floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dVelocityGradientdStateVariables,                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousVelocityGradient,                                    floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousCauchyStress,              floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousF,                         floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousSubFs,                     floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousT,                         floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousVelocityGradientdPreviousStateVariables,            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableEvolutionRates,                                 floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdCauchyStress,                   floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdF,                              floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdSubFs,                          floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdT,                              floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableEvolutionRatesdStateVariables,                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariableEvolutionRates,                         floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousCauchyStress,   floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousF,              floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousSubFs,          floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousT,              floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, dPreviousStateVariableEvolutionRatesdPreviousStateVariables, floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticDeformationGradient,                                  floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdCauchyStress,                    floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdF,                               floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdSubFs,                           floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdT,                               floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdStateVariables,                  floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousCauchyStress,            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousF,                       floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousSubFs,                   floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousT,                       floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticDeformationGradientdPreviousStateVariables,          floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStateVariables,                                       floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdCauchyStress,                         floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdF,                                    floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdSubFs,                                floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdT,                                    floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdStateVariables,                       floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousCauchyStress,                 floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousF,                            floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousSubFs,                        floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousT,                            floatVector )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStateVariablesdPreviousStateVariables,               floatMatrix )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariables,                                              floatVector )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousStateVariables,                                      floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, peryznaParameters,                                           floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, dragStressParameters,                                        floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, thermalParameters,                                           floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, yieldParameters,                                             floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, flowParameters,                                              floatVector )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, hardeningParameters,                                         floatVector )

        };

    }

}

#endif
