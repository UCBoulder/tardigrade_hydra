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

                void setFlowDirection( const floatVector &flowDirection );

                void setPreviousFlowDirection( const floatVector &previousFlowDirection );

                void setYieldFunction( const floatType &yieldFunction );

                void setPreviousYieldFunction( const floatType &previousYieldFunction );

                void setPlasticThermalMultiplier( const floatType &plasticThermalMultiplier );

                void setPreviousPlasticThermalMultiplier( const floatType &previousPlasticThermalMultiplier );

                void setDragStress( const floatType &dragStress );

                void setPreviousDragStress( const floatType &previousDragStress );

                void setHardeningFunction( const floatType &hardeningFunction );

                void setPreviousHardeningFunction( const floatType &previousHardeningFunction );

                void setPlasticMultiplier( const floatType &plasticMultiplier );

                void setPreviousPlasticMultiplier( const floatType &previousPlasticMultiplier );

                void setVelocityGradient( const floatVector &velocityGradient );

                void setPreviousVelocityGradient( const floatVector &previousVelocityGradient );

                void setStateVariableEvolutionRates( const floatVector &stateVariableEvolutionRate );

                void setPreviousStateVariableEvolutionRates( const floatVector &previousStateVariableEvolutionRates );

                void setPlasticDeformationGradient( const floatVector &plasticDeformationGradient );

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

                const floatVector* getFlowDirection( );

                const floatType* getYieldFunction( );

                const floatType* getPlasticThermalMultiplier( );

                const floatType* getDragStress( );

                const floatType* getHardeningFunction( );

                const floatType* getPlasticMultiplier( );

                const floatVector* getVelocityGradient( );

                const floatVector* getStateVariableEvolutionRates( );

                const floatVector* getPreviousDrivingStress( );

                const floatVector* getPreviousFlowDirection( );

                const floatType* getPreviousYieldFunction( );

                const floatType* getPreviousPlasticThermalMultiplier( );

                const floatType* getPreviousDragStress( );

                const floatType* getPreviousHardeningFunction( );

                const floatType* getPreviousPlasticMultiplier( );

                const floatVector* getPreviousVelocityGradient( );

                const floatVector* getPreviousStateVariableEvolutionRates( );

                const floatVector* getPlasticDeformationGradient( );

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

                virtual void setFlowDirection( );

                virtual void setYieldFunction( );

                virtual void setPlasticThermalMultiplier( );

                virtual void setDragStress( );

                virtual void setHardeningFunction( );

                virtual void setPlasticMultiplier( );

                virtual void setVelocityGradient( );

                virtual void setStateVariableEvolutionRates( );

                virtual void setPreviousDrivingStress( );

                virtual void setPreviousFlowDirection( );

                virtual void setPreviousYieldFunction( );

                virtual void setPreviousPlasticThermalMultiplier( );

                virtual void setPreviousDragStress( );

                virtual void setPreviousHardeningFunction( );

                virtual void setPreviousPlasticMultiplier( );

                virtual void setPreviousVelocityGradient( );

                virtual void setPreviousStateVariableEvolutionRates( );

                virtual void setDrivingStress( const bool isPrevious );

                virtual void setFlowDirection( const bool isPrevious );

                virtual void setYieldFunction( const bool isPrevious);

                virtual void setPlasticThermalMultiplier( const bool isPrevious );

                virtual void setDragStress( const bool isPrevious );

                virtual void setHardeningFunction( const bool isPrevious );

                virtual void setPlasticMultiplier( const bool isPrevious );

                virtual void setVelocityGradient( const bool isPrevious );

                virtual void setStateVariableEvolutionRates( const bool isPrevious );

                virtual void setPlasticDeformationGradient( );

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

                tardigradeHydra::dataStorage< floatVector > _flowDirection;

                tardigradeHydra::dataStorage< floatVector > _previousFlowDirection;

                tardigradeHydra::dataStorage< floatType > _yieldFunction;

                tardigradeHydra::dataStorage< floatType > _previousYieldFunction;

                tardigradeHydra::dataStorage< floatType > _plasticThermalMultiplier;

                tardigradeHydra::dataStorage< floatType > _previousPlasticThermalMultiplier;

                tardigradeHydra::dataStorage< floatType > _dragStress;

                tardigradeHydra::dataStorage< floatType > _previousDragStress;

                tardigradeHydra::dataStorage< floatType > _hardeningFunction;

                tardigradeHydra::dataStorage< floatType > _previousHardeningFunction;

                tardigradeHydra::dataStorage< floatType > _plasticMultiplier;

                tardigradeHydra::dataStorage< floatType > _previousPlasticMultiplier;

                tardigradeHydra::dataStorage< floatVector > _velocityGradient;

                tardigradeHydra::dataStorage< floatVector > _previousVelocityGradient;

                tardigradeHydra::dataStorage< floatVector > _stateVariableEvolutionRates;

                tardigradeHydra::dataStorage< floatVector > _previousStateVariableEvolutionRates;

                tardigradeHydra::dataStorage< floatVector > _plasticDeformationGradient;

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
