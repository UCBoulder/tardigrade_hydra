/**
  ******************************************************************************
  * \file tardigrade-hydraPeryznaViscoplasticity.h
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

                void setDrivingStress( const floatVector &drivingStress );

                void setPreviousDrivingStress( const floatVector &previousDrivingStress );

                void setdDrivingStressdCauchyStress( const floatMatrix &dDrivingStressdCauchyStress );

                void setdDrivingStressdF( const floatMatrix &dDrivingStressdF );

                void setdDrivingStressdSubFs( const floatMatrix &dDrivingStressdSubFs );

                void setdPreviousDrivingStressdPreviousCauchyStress( const floatMatrix &dPreviousDrivingStressdPreviousCauchyStress );

                void setdPreviousDrivingStressdPreviousF( const floatMatrix &dPreviousDrivingStressdPreviousF );

                void setdPreviousDrivingStressdPreviousSubFs( const floatMatrix &dPreviousDrivingStressdPreviousSubFs );

                void setFlowDirection( const floatVector &flowDirection );

                void setPreviousFlowDirection( const floatVector &previousFlowDirection );

                void setdFlowDirectiondCauchyStress( const floatMatrix &dFlowDirectiondCauchyStress );

                void setdFlowDirectiondF( const floatMatrix &dFlowDirectiondF );

                void setdFlowDirectiondSubFs( const floatMatrix &dFlowDirectiondSubFs );

                void setdPreviousFlowDirectiondPreviousCauchyStress( const floatMatrix &dPreviousFlowDirectiondPreviousCauchyStress );

                void setdPreviousFlowDirectiondPreviousF( const floatMatrix &dPreviousFlowDirectiondPreviousF );

                void setdPreviousFlowDirectiondPreviousSubFs( const floatMatrix &dPreviousFlowDirectiondPreviousSubFs );

                void setYieldFunction( const floatType &yieldFunction );

                void setPreviousYieldFunction( const floatType &previousYieldFunction );

                void setdYieldFunctiondCauchyStress( const floatVector &dYieldFunctiondCauchyStress );

                void setdYieldFunctiondF( const floatVector &dYieldFunctiondF );

                void setdYieldFunctiondSubFs( const floatVector &dYieldFunctiondSubFs );

                void setdPreviousYieldFunctiondPreviousCauchyStress( const floatVector &dPreviousYieldFunctiondPreviousCauchyStress );

                void setdPreviousYieldFunctiondPreviousF( const floatVector &dPreviousYieldFunctiondPreviousF );

                void setdPreviousYieldFunctiondPreviousSubFs( const floatVector &dPreviousYieldFunctiondPreviousSubFs );

                void setPlasticThermalMultiplier( const floatType &plasticThermalMultiplier );

                void setPreviousPlasticThermalMultiplier( const floatType &previousPlasticThermalMultiplier );

                void setdPlasticThermalMultiplierdT( const floatType &dPlasticThermalMultiplierdT );

                void setdPreviousPlasticThermalMultiplierdPreviousT( const floatType &dPreviousPlasticThermalMultiplierdPreviousT );

                void setDragStress( const floatType &dragStress );

                void setPreviousDragStress( const floatType &previousDragStress );

                void setdDragStressdStateVariables( const floatVector &dDragStressdStateVariables );

                void setdPreviousDragStressdPreviousStateVariables( const floatVector &dPreviousDragStressdPreviousStateVariables );

                void setHardeningFunction( const floatType &hardeningFunction );

                void setdHardeningFunctiondStateVariables( const floatVector &dHardeningFunctiondStateVariables );

                void setPreviousHardeningFunction( const floatType &previousHardeningFunction );

                void setdPreviousHardeningFunctiondPreviousStateVariables( const floatVector &dPreviousHardeningFunctiondPreviousStateVariables );

                void setPlasticMultiplier( const floatType &plasticMultiplier );

                void setdPlasticMultiplierdCauchyStress( const floatVector &dPlasticMultiplierdCauchyStress );

                void setdPlasticMultiplierdF( const floatVector &dPlasticMultiplierdF );

                void setdPlasticMultiplierdSubFs( const floatVector &dPlasticMultiplierdSubFs );

                void setdPlasticMultiplierdT( const floatType &dPlasticMultiplierdT );

                void setdPlasticMultiplierdStateVariables( const floatVector &dPlasticMultiplierdStateVariables );

                void setPreviousPlasticMultiplier( const floatType &previousPlasticMultiplier );

                void setdPreviousPlasticMultiplierdPreviousCauchyStress( const floatVector &dPreviousPlasticMultiplierdPreviousCauchyStress );

                void setdPreviousPlasticMultiplierdPreviousF( const floatVector &dPreviousPlasticMultiplierdPreviousF );

                void setdPreviousPlasticMultiplierdPreviousSubFs( const floatVector &dPreviousPlasticMultiplierdPreviousSubFs );

                void setdPreviousPlasticMultiplierdPreviousT( const floatType &dPreviousPlasticMultiplierdPreviousT );

                void setdPreviousPlasticMultiplierdPreviousStateVariables( const floatVector &dPreviousPlasticMultiplierdPreviousStateVariables );

                void setVelocityGradient( const floatVector &velocityGradient );

                void setdVelocityGradientdCauchyStress( const floatMatrix &dVelocityGradientdCauchyStress );

                void setdVelocityGradientdF( const floatMatrix &dVelocityGradientdF );

                void setdVelocityGradientdSubFs( const floatMatrix &dVelocityGradientdSubFs );

                void setdVelocityGradientdT( const floatVector &dVelocityGradientdT );

                void setdVelocityGradientdStateVariables( const floatMatrix &dVelocityGradientdStateVariables );

                void setPreviousVelocityGradient( const floatVector &previousVelocityGradient );

                void setdPreviousVelocityGradientdPreviousCauchyStress( const floatMatrix &dPreviousVelocityGradientdPreviousCauchyStress );

                void setdPreviousVelocityGradientdPreviousF( const floatMatrix &dPreviousVelocityGradientdPreviousF );

                void setdPreviousVelocityGradientdPreviousSubFs( const floatMatrix &dPreviousVelocityGradientdPreviousSubFs );

                void setdPreviousVelocityGradientdPreviousT( const floatVector &dPreviousVelocityGradientdPreviousT );

                void setdPreviousVelocityGradientdPreviousStateVariables( const floatMatrix &dPreviousVelocityGradientdPreviousStateVariables );

                void setStateVariableEvolutionRates( const floatVector &stateVariableEvolutionRates );

                void setdStateVariableEvolutionRatesdCauchyStress( const floatMatrix &dStateVariableEvolutionRatesdCauchyStress );

                void setdStateVariableEvolutionRatesdF( const floatMatrix &dStateVariableEvolutionRatesdF );

                void setdStateVariableEvolutionRatesdSubFs( const floatMatrix &dStateVariableEvolutionRatesdSubFs );

                void setdStateVariableEvolutionRatesdT( const floatVector &dStateVariableEvolutionRatesdT );

                void setdStateVariableEvolutionRatesdStateVariables( const floatMatrix &dStateVariableEvolutionRatesdStateVariables );

                void setPreviousStateVariableEvolutionRates( const floatVector &previousStateVariableEvolutionRates );

                void setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousCauchyStress );

                void setdPreviousStateVariableEvolutionRatesdPreviousF( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousF );

                void setdPreviousStateVariableEvolutionRatesdPreviousSubFs( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousSubFs );

                void setdPreviousStateVariableEvolutionRatesdPreviousT( const floatVector &dPreviousStateVariableEvolutionRatesdPreviousT );

                void setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousStateVariables );

                void setPlasticDeformationGradient( const floatVector &plasticDeformationGradient );

                void setdPlasticDeformationGradientdCauchyStress( const floatMatrix &dPlasticDeformationGradientdCauchyStress );

                void setdPlasticDeformationGradientdF( const floatMatrix &dPlasticDeformationGradientdF );

                void setdPlasticDeformationGradientdSubFs( const floatMatrix &dPlasticDeformationGradientdSubFs );

                void setdPlasticDeformationGradientdT( const floatVector &dPlasticDeformationGradientdT );

                void setdPlasticDeformationGradientdStateVariables( const floatMatrix &dPlasticDeformationGradientdStateVariables );

                void setdPlasticDeformationGradientdPreviousCauchyStress( const floatMatrix &dPlasticDeformationGradientdPreviousCauchyStress );

                void setdPlasticDeformationGradientdPreviousF( const floatMatrix &dPlasticDeformationGradientdPreviousF );

                void setdPlasticDeformationGradientdPreviousSubFs( const floatMatrix &dPlasticDeformationGradientdPreviousSubFs );

                void setdPlasticDeformationGradientdPreviousT( const floatVector &dPlasticDeformationGradientdPreviousT );

                void setdPlasticDeformationGradientdPreviousStateVariables( const floatMatrix &dPlasticDeformationGradientdPreviousStateVariables );

                void setPlasticStateVariables( const floatVector &plasticStateVariables );

                void setStateVariables( const floatVector &stateVariables );

                void setPreviousStateVariables( const floatVector &previousStateVariables );

                void setPeryznaParameters( const floatVector &peryznaParameters );

                void setDragStressParameters( const floatVector &dragStressParameters );

                void setThermalParameters( const floatVector &thermalParameters );

                void setYieldParameters( const floatVector &yieldParameters );

                void setFlowParameters( const floatVector &flowParameters );

                void setHardeningParameters( const floatVector &hardeningParameters );

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

                const floatVector* getStateVariables( );

                const floatVector* getPreviousStateVariables( );

                const floatVector* getPeryznaParameters( );

                const floatVector* getDragStressParameters( );

                const floatVector* getThermalParameters( );

                const floatVector* getYieldParameters( );

                const floatVector* getFlowParameters( );

                const floatVector* getHardeningParameters( );

                const floatType* getIntegrationParameter( );

            private:

                unsigned int _plasticConfigurationIndex;

                std::vector< unsigned int > _stateVariableIndices;

                floatType _integrationParameter;

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

                virtual void setStateVariables( );

                virtual void setStateVariables( const bool isPrevious );

                virtual void setPreviousStateVariables( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdF( ) override;

                virtual void setdRdT( ) override;

                virtual void decomposeParameters( const floatVector &parameters );

                tardigradeHydra::dataStorage< floatVector > _drivingStress;

                tardigradeHydra::dataStorage< floatVector > _previousDrivingStress;

                tardigradeHydra::dataStorage< floatMatrix > _dDrivingStressdCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dDrivingStressdF;

                tardigradeHydra::dataStorage< floatMatrix > _dDrivingStressdSubFs;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousDrivingStressdPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousDrivingStressdPreviousF;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousDrivingStressdPreviousSubFs;

                tardigradeHydra::dataStorage< floatVector > _flowDirection;

                tardigradeHydra::dataStorage< floatVector > _previousFlowDirection;

                tardigradeHydra::dataStorage< floatMatrix > _dFlowDirectiondCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dFlowDirectiondF;

                tardigradeHydra::dataStorage< floatMatrix > _dFlowDirectiondSubFs;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousFlowDirectiondPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousFlowDirectiondPreviousF;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousFlowDirectiondPreviousSubFs;

                tardigradeHydra::dataStorage< floatType > _yieldFunction;

                tardigradeHydra::dataStorage< floatVector > _dYieldFunctiondCauchyStress;

                tardigradeHydra::dataStorage< floatVector > _dYieldFunctiondF;

                tardigradeHydra::dataStorage< floatVector > _dYieldFunctiondSubFs;

                tardigradeHydra::dataStorage< floatType > _previousYieldFunction;

                tardigradeHydra::dataStorage< floatVector > _dPreviousYieldFunctiondPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatVector > _dPreviousYieldFunctiondPreviousF;

                tardigradeHydra::dataStorage< floatVector > _dPreviousYieldFunctiondPreviousSubFs;

                tardigradeHydra::dataStorage< floatType > _plasticThermalMultiplier;

                tardigradeHydra::dataStorage< floatType > _dPlasticThermalMultiplierdT;

                tardigradeHydra::dataStorage< floatType > _previousPlasticThermalMultiplier;

                tardigradeHydra::dataStorage< floatType > _dPreviousPlasticThermalMultiplierdPreviousT;

                tardigradeHydra::dataStorage< floatType > _dragStress;

                tardigradeHydra::dataStorage< floatVector > _dDragStressdStateVariables;

                tardigradeHydra::dataStorage< floatType > _previousDragStress;

                tardigradeHydra::dataStorage< floatVector > _dPreviousDragStressdPreviousStateVariables;

                tardigradeHydra::dataStorage< floatType > _hardeningFunction;

                tardigradeHydra::dataStorage< floatVector > _dHardeningFunctiondStateVariables;

                tardigradeHydra::dataStorage< floatType > _previousHardeningFunction;

                tardigradeHydra::dataStorage< floatVector > _dPreviousHardeningFunctiondPreviousStateVariables;

                tardigradeHydra::dataStorage< floatType > _plasticMultiplier;

                tardigradeHydra::dataStorage< floatVector > _dPlasticMultiplierdCauchyStress;

                tardigradeHydra::dataStorage< floatVector > _dPlasticMultiplierdF;

                tardigradeHydra::dataStorage< floatVector > _dPlasticMultiplierdSubFs;

                tardigradeHydra::dataStorage< floatType > _dPlasticMultiplierdT;

                tardigradeHydra::dataStorage< floatVector > _dPlasticMultiplierdStateVariables;

                tardigradeHydra::dataStorage< floatType > _previousPlasticMultiplier;

                tardigradeHydra::dataStorage< floatVector > _dPreviousPlasticMultiplierdPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatVector > _dPreviousPlasticMultiplierdPreviousF;

                tardigradeHydra::dataStorage< floatVector > _dPreviousPlasticMultiplierdPreviousSubFs;

                tardigradeHydra::dataStorage< floatType > _dPreviousPlasticMultiplierdPreviousT;

                tardigradeHydra::dataStorage< floatVector > _dPreviousPlasticMultiplierdPreviousStateVariables;

                tardigradeHydra::dataStorage< floatVector > _velocityGradient;

                tardigradeHydra::dataStorage< floatMatrix > _dVelocityGradientdCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dVelocityGradientdF;

                tardigradeHydra::dataStorage< floatMatrix > _dVelocityGradientdSubFs;

                tardigradeHydra::dataStorage< floatVector > _dVelocityGradientdT;

                tardigradeHydra::dataStorage< floatMatrix > _dVelocityGradientdStateVariables;

                tardigradeHydra::dataStorage< floatVector > _previousVelocityGradient;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousVelocityGradientdPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousVelocityGradientdPreviousF;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousVelocityGradientdPreviousSubFs;

                tardigradeHydra::dataStorage< floatVector > _dPreviousVelocityGradientdPreviousT;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousVelocityGradientdPreviousStateVariables;

                tardigradeHydra::dataStorage< floatVector > _stateVariableEvolutionRates;

                tardigradeHydra::dataStorage< floatMatrix > _dStateVariableEvolutionRatesdCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dStateVariableEvolutionRatesdF;

                tardigradeHydra::dataStorage< floatMatrix > _dStateVariableEvolutionRatesdSubFs;

                tardigradeHydra::dataStorage< floatVector > _dStateVariableEvolutionRatesdT;

                tardigradeHydra::dataStorage< floatMatrix > _dStateVariableEvolutionRatesdStateVariables;

                tardigradeHydra::dataStorage< floatVector > _previousStateVariableEvolutionRates;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousStateVariableEvolutionRatesdPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousStateVariableEvolutionRatesdPreviousF;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousStateVariableEvolutionRatesdPreviousSubFs;

                tardigradeHydra::dataStorage< floatVector > _dPreviousStateVariableEvolutionRatesdPreviousT;

                tardigradeHydra::dataStorage< floatMatrix > _dPreviousStateVariableEvolutionRatesdPreviousStateVariables;

                tardigradeHydra::dataStorage< floatVector > _plasticDeformationGradient;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdF;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdSubFs;

                tardigradeHydra::dataStorage< floatVector > _dPlasticDeformationGradientdT;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdStateVariables;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdPreviousCauchyStress;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdPreviousF;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdPreviousSubFs;

                tardigradeHydra::dataStorage< floatVector > _dPlasticDeformationGradientdPreviousT;

                tardigradeHydra::dataStorage< floatMatrix > _dPlasticDeformationGradientdPreviousStateVariables;

                tardigradeHydra::dataStorage< floatVector > _plasticStateVariables;

                tardigradeHydra::dataStorage< floatVector > _stateVariables;

                tardigradeHydra::dataStorage< floatVector > _previousStateVariables;

                tardigradeHydra::dataStorage< floatVector > _peryznaParameters;

                tardigradeHydra::dataStorage< floatVector > _dragStressParameters;

                tardigradeHydra::dataStorage< floatVector > _thermalParameters;

                tardigradeHydra::dataStorage< floatVector > _yieldParameters;

                tardigradeHydra::dataStorage< floatVector > _flowParameters;

                tardigradeHydra::dataStorage< floatVector > _hardeningParameters;

        };

    }

}

#endif
