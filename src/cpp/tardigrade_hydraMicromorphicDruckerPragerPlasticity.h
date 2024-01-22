/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicLinearElasticity.h
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework. Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_MICROMORPHIC_DRUCKER_PRAGER_PLASTICITY_H
#define TARDIGRADE_HYDRA_MICROMORPHIC_DRUCKER_PRAGER_PLASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydraMicromorphic.h>

namespace tardigradeHydra{

    namespace micromorphicDruckerPragerPlasticity{

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

        typedef floatType variableType; //!< Define the variable values type.
        typedef std::vector< variableType > variableVector; //!< Define a vector of variables
        typedef std::vector< std::vector< variableType > > variableMatrix; //!< Define a matrix of variables

        typedef double parameterType; //!< Define the parameter values type.
        typedef std::vector< parameterType > parameterVector; //!< Define a vector of parameters

        typedef double constantType; //!< Define the constant values type.
        typedef std::vector< constantType > constantVector; //!< Define a vector of constants
        typedef std::vector< std::vector< constantType > > constantMatrix; //!< Define a matrix of constants

        void computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B );

        void computeSecondOrderDruckerPragerYieldEquation( const floatVector   &stressMeasure,
                                                           const floatType     &cohesion,
                                                           const floatVector   &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue );

        void computeSecondOrderDruckerPragerYieldEquation( const floatVector   &stressMeasure,
                                                           const floatType     &cohesion,
                                                           const floatVector   &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue, floatVector &dFdStress, floatType &dFdc,
                                                           floatVector &dFdPrecedingF, double tol = 1e-9 );

        void computeSecondOrderDruckerPragerYieldEquation( const floatVector   &stressMeasure,
                                                           const floatType     &cohesion,
                                                           const floatVector   &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue, floatVector &dFdStress, floatType &dFdc,
                                                           floatVector &dFdPrecedingF, floatMatrix &d2FdStress2,
                                                           floatMatrix &d2FdStressdPrecedingF, double tol = 1e-9 );

        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                           const variableVector &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue );
    
        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                           const variableVector &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdPrecedingF );
    
        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                           const variableVector &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdPrecedingF, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdPrecedingF );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma,
                                                  variableMatrix &dPlasticMacroLdElasticRCG,
                                                  variableMatrix &dPlasticMacroLdMacroFlowDirection,
                                                  variableMatrix &dPlasticMacroLdMicroFlowDirection );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma,
                                                  variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                                  variableMatrix &dPlasticMicroLdElasticPsi,
                                                  variableMatrix &dPlasticMicroLdMicroFlowDirection );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                                          variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableVector &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                                          variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &LHS,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                        variableMatrix &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableMatrix &dPlasticFdPlasticMacroL,
                                       variableMatrix &dPlasticMicroDeformationdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMacroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableMatrix &dPlasticFdPlasticMacroL,
                                       variableMatrix &dPlasticMicroDeformationdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMacroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL,
                                       variableMatrix &dPlasticFdPreviousPlasticF,
                                       variableMatrix &dPlasticFdPreviousPlasticMacroL,
                                       variableMatrix &dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                       variableMatrix &dPlasticMicroDeformationdPreviousPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                       variableMatrix &dPlasticMicroGradientdPreviousPlasticMicroGradient,
                                       variableMatrix &dPlasticMicroGradientdPreviousPlasticMacroL,
                                       variableMatrix &dPlasticMicroGradientdPreviousPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPreviousPlasticMicroGradientL,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        /*!
         * The residual for a micromorphic Drucker Prager plasticity model
         */
        class residual : public tardigradeHydra::residualBaseMicromorphic {

            public:

                using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

                using tardigradeHydra::residualBaseMicromorphic::setResidual;

                using tardigradeHydra::residualBaseMicromorphic::setJacobian;

                using tardigradeHydra::residualBaseMicromorphic::setdRdD;

                using tardigradeHydra::residualBaseMicromorphic::setdRdT;

                residual( hydraBaseMicromorphic *_hydra, const unsigned int &_numEquations, const unsigned int &plasticConfigurationIndex,
                          const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const floatType integrationParameter = 0.5 )
                        : tardigradeHydra::residualBaseMicromorphic( _hydra, _numEquations ){
                    /*!
                     * The main initialization constructor for the Drucker Prager plasticity residual
                     * 
                     * \param *_hydra: A pointer to the containing hydra class
                     * \param &_numEquations: The number of equations the residual defines
                     * \param &plasticConfigurationIndex: The index of the configuration which represents the plastic deformation
                     * \param &stateVariableIndices: The indices of the plastic state variables
                     * \param &parameters: The parameter vector
                     * \param &integrationParameter: The integration parameter for the function. 0 is explicit, 1 is implicit.
                     */

                    _plasticConfigurationIndex = plasticConfigurationIndex;

                    _numPlasticMultipliers = 5; //Hard-coding this because we couldn't handle a different number anyway.

                    _stateVariableIndices = stateVariableIndices;

                    _integrationParameter = integrationParameter;

                    TARDIGRADE_ERROR_TOOLS_CATCH( extractMaterialParameters( parameters ) );

                }

                //!Get the number of plastic multipliers expected for the problem
                const unsigned int* getNumPlasticMultipliers( ){ return &_numPlasticMultipliers; }

                const unsigned int* getPlasticConfigurationIndex( );

                const std::vector< unsigned int >* getStateVariableIndices( );

                const floatType* getIntegrationParameter( );

            protected:

                unsigned int _numPlasticMultipliers; //!< The number of plastic multipliers. Hard coded to 5 but setting as a variable just in case

                virtual void extractMaterialParameters( const parameterVector &parameters );

                virtual void setMacroDrivingStress( );

                virtual void setSymmetricMicroDrivingStress( );

                virtual void setHigherOrderDrivingStress( );

                virtual void setPreviousMacroDrivingStress( );

                virtual void setPreviousSymmetricMicroDrivingStress( );

                virtual void setPreviousHigherOrderDrivingStress( );

                virtual void setDrivingStresses( const bool isPrevious );

                virtual void setdMacroDrivingStressdMacroStress( );

                virtual void setdSymmetricMicroDrivingStressdMicroStress( );

                virtual void setdHigherOrderDrivingStressdHigherOrderStress( );

                virtual void setdMacroDrivingStressdF( );

                virtual void setdSymmetricMicroDrivingStressdF( );

                virtual void setdHigherOrderDrivingStressdF( );

                virtual void setdHigherOrderDrivingStressdChi( );

                virtual void setdMacroDrivingStressdFn( );

                virtual void setdSymmetricMicroDrivingStressdFn( );

                virtual void setdHigherOrderDrivingStressdFn( );

                virtual void setdHigherOrderDrivingStressdChin( );

                virtual void setPreviousdMacroDrivingStressdMacroStress( );

                virtual void setPreviousdSymmetricMicroDrivingStressdMicroStress( );

                virtual void setPreviousdHigherOrderDrivingStressdHigherOrderStress( );

                virtual void setPreviousdMacroDrivingStressdF( );

                virtual void setPreviousdSymmetricMicroDrivingStressdF( );

                virtual void setPreviousdHigherOrderDrivingStressdF( );

                virtual void setPreviousdHigherOrderDrivingStressdChi( );

                virtual void setPreviousdMacroDrivingStressdFn( );

                virtual void setPreviousdSymmetricMicroDrivingStressdFn( );

                virtual void setPreviousdHigherOrderDrivingStressdFn( );

                virtual void setPreviousdHigherOrderDrivingStressdChin( );

                virtual void setDrivingStressesJacobians( const bool isPrevious );

                virtual void setPlasticStateVariables( );

                virtual void setPreviousPlasticStateVariables( );

                virtual void setPlasticStateVariables( const bool isPrevious );

                virtual void setPlasticMultipliers( );

                virtual void setPreviousPlasticMultipliers( );

                virtual void setPlasticMultipliers( const bool isPrevious );

                virtual void setPlasticStrainLikeISVs( );

                virtual void setPreviousPlasticStrainLikeISVs( );

                virtual void setPlasticStrainLikeISVs( const bool isPrevious );

                virtual void setMacroCohesion( );

                virtual void setMicroCohesion( );

                virtual void setMicroGradientCohesion( );

                virtual void setPreviousMacroCohesion( );

                virtual void setPreviousMicroCohesion( );

                virtual void setPreviousMicroGradientCohesion( );

                virtual void setCohesions( const bool isPrevious );

                virtual void setdMacroCohesiondStateVariables( );

                virtual void setdMicroCohesiondStateVariables( );

                virtual void setdMicroGradientCohesiondStateVariables( );

                virtual void setPreviousdMacroCohesiondStateVariables( );

                virtual void setPreviousdMicroCohesiondStateVariables( );

                virtual void setPreviousdMicroGradientCohesiondStateVariables( );

                virtual void setCohesionsJacobians( const bool isPrevious );

                virtual void setdMacroFlowdc( );

                virtual void setdMicroFlowdc( );

                virtual void setdMicroGradientFlowdc( );

                virtual void setdMacroFlowdDrivingStress( );

                virtual void setdMicroFlowdDrivingStress( );

                virtual void setdMicroGradientFlowdDrivingStress( );

                virtual void setPreviousdMacroFlowdc( );

                virtual void setPreviousdMicroFlowdc( );

                virtual void setPreviousdMicroGradientFlowdc( );

                virtual void setPreviousdMacroFlowdDrivingStress( );

                virtual void setPreviousdMicroFlowdDrivingStress( );

                virtual void setPreviousdMicroGradientFlowdDrivingStress( );

                virtual void setFlowPotentialGradients( const bool isPrevious );

                virtual void setd2MacroFlowdDrivingStressdStress( );

                virtual void setd2MacroFlowdDrivingStressdF( );

                virtual void setd2MacroFlowdDrivingStressdFn( );

                virtual void setd2MicroFlowdDrivingStressdStress( );

                virtual void setd2MicroFlowdDrivingStressdF( );

                virtual void setd2MicroFlowdDrivingStressdFn( );

                virtual void setd2MicroGradientFlowdDrivingStressdStress( );

                virtual void setd2MicroGradientFlowdDrivingStressdF( );

                virtual void setd2MicroGradientFlowdDrivingStressdFn( );

                virtual void setd2MicroGradientFlowdDrivingStressdChi( );

                virtual void setd2MicroGradientFlowdDrivingStressdChin( );

                virtual void setPreviousd2MacroFlowdDrivingStressdStress( );

                virtual void setPreviousd2MacroFlowdDrivingStressdF( );

                virtual void setPreviousd2MacroFlowdDrivingStressdFn( );

                virtual void setPreviousd2MicroFlowdDrivingStressdStress( );

                virtual void setPreviousd2MicroFlowdDrivingStressdF( );

                virtual void setPreviousd2MicroFlowdDrivingStressdFn( );

                virtual void setPreviousd2MicroGradientFlowdDrivingStressdStress( );

                virtual void setPreviousd2MicroGradientFlowdDrivingStressdF( );

                virtual void setPreviousd2MicroGradientFlowdDrivingStressdFn( );

                virtual void setPreviousd2MicroGradientFlowdDrivingStressdChi( );

                virtual void setPreviousd2MicroGradientFlowdDrivingStressdChin( );

                virtual void setFlowPotentialGradientsJacobians( const bool isPrevious );

                virtual void setPlasticStrainLikeISVEvolutionRates( );

                virtual void setPreviousPlasticStrainLikeISVEvolutionRates( );

                virtual void setPlasticStrainLikeISVEvolutionRates( const bool isPrevious );

                virtual void setdPlasticStrainLikeISVEvolutionRatesdStateVariables( );

                virtual void setPreviousdPlasticStrainLikeISVEvolutionRatesdStateVariables( );

                virtual void setPlasticStrainLikeISVEvolutionRatesJacobians( const bool isPrevious );

                virtual void setUpdatedPlasticStrainLikeISVs( );

                virtual void setdUpdatedPlasticStrainLikeISVsdStateVariables( );

                virtual void setdUpdatedPlasticStrainLikeISVsdPreviousStateVariables( );

                virtual void setUpdatedPlasticStrainLikeISVsJacobians( const bool addPrevious );

                virtual void setMacroYield( );

                virtual void setMicroYield( );

                virtual void setMicroGradientYield( );

                virtual void setPreviousMacroYield( );

                virtual void setPreviousMicroYield( );

                virtual void setPreviousMicroGradientYield( );

                virtual void setYield( const bool isPrevious );

                virtual void setdMacroYielddStress( );

                virtual void setdMacroYielddStateVariables( );

                virtual void setdMacroYielddF( );

                virtual void setdMacroYielddFn( );

                virtual void setdMicroYielddStress( );

                virtual void setdMicroYielddStateVariables( );

                virtual void setdMicroYielddF( );

                virtual void setdMicroYielddFn( );

                virtual void setdMicroGradientYielddStress( );

                virtual void setdMicroGradientYielddStateVariables( );

                virtual void setdMicroGradientYielddF( );

                virtual void setdMicroGradientYielddFn( );

                virtual void setdMicroGradientYielddChi( );

                virtual void setdMicroGradientYielddChin( );

                virtual void setPreviousdMacroYielddStress( );

                virtual void setPreviousdMacroYielddStateVariables( );

                virtual void setPreviousdMacroYielddF( );

                virtual void setPreviousdMacroYielddFn( );

                virtual void setPreviousdMicroYielddStress( );

                virtual void setPreviousdMicroYielddStateVariables( );

                virtual void setPreviousdMicroYielddF( );

                virtual void setPreviousdMicroYielddFn( );

                virtual void setPreviousdMicroGradientYielddStress( );

                virtual void setPreviousdMicroGradientYielddStateVariables( );

                virtual void setPreviousdMicroGradientYielddF( );

                virtual void setPreviousdMicroGradientYielddFn( );

                virtual void setPreviousdMicroGradientYielddChi( );

                virtual void setPreviousdMicroGradientYielddChin( );

                virtual void setYieldJacobians( const bool isPrevious );

                virtual void setPrecedingDeformationGradient( );

                virtual void setPreviousPrecedingDeformationGradient( );

                virtual void setPrecedingDeformationGradient( const bool isPrevious );

                virtual void setdPrecedingDeformationGradientdF( );

                virtual void setdPrecedingDeformationGradientdFn( );

                virtual void setPreviousdPrecedingDeformationGradientdF( );

                virtual void setPreviousdPrecedingDeformationGradientdFn( );

                virtual void setPrecedingDeformationGradientJacobians( const bool isPrevious );

                virtual void setPrecedingMicroDeformation( );

                virtual void setPreviousPrecedingMicroDeformation( );

                virtual void setPrecedingMicroDeformation( const bool isPrevious );

                virtual void setdPrecedingMicroDeformationdChi( );

                virtual void setdPrecedingMicroDeformationdChin( );

                virtual void setPreviousdPrecedingMicroDeformationdChi( );

                virtual void setPreviousdPrecedingMicroDeformationdChin( );

                virtual void setPrecedingMicroDeformationJacobians( const bool isPrevious );

                virtual void setPlasticMacroVelocityGradient( );

                virtual void setPreviousPlasticMacroVelocityGradient( );

                virtual void setPlasticMicroVelocityGradient( );

                virtual void setPreviousPlasticMicroVelocityGradient( );

                virtual void setPlasticGradientMicroVelocityGradient( );

                virtual void setPreviousPlasticGradientMicroVelocityGradient( );

                virtual void setPlasticVelocityGradients( const bool isPrevious );

                virtual void setdPlasticMacroVelocityGradientdMacroStress( );

                virtual void setPreviousdPlasticMacroVelocityGradientdMacroStress( );

                virtual void setdPlasticMacroVelocityGradientdMicroStress( );

                virtual void setPreviousdPlasticMacroVelocityGradientdMicroStress( );

                virtual void setdPlasticMacroVelocityGradientdF( );

                virtual void setPreviousdPlasticMacroVelocityGradientdF( );

                virtual void setdPlasticMacroVelocityGradientdFn( );

                virtual void setPreviousdPlasticMacroVelocityGradientdFn( );

                virtual void setdPlasticMacroVelocityGradientdStateVariables( );

                virtual void setPreviousdPlasticMacroVelocityGradientdStateVariables( );

                virtual void setdPlasticMicroVelocityGradientdMicroStress( );

                virtual void setPreviousdPlasticMicroVelocityGradientdMicroStress( );

                virtual void setdPlasticMicroVelocityGradientdF( );

                virtual void setPreviousdPlasticMicroVelocityGradientdF( );

                virtual void setdPlasticMicroVelocityGradientdFn( );

                virtual void setPreviousdPlasticMicroVelocityGradientdFn( );

                virtual void setdPlasticMicroVelocityGradientdChi( );

                virtual void setdPlasticMicroVelocityGradientdChin( );

                virtual void setPreviousdPlasticMicroVelocityGradientdChi( );

                virtual void setPreviousdPlasticMicroVelocityGradientdChin( );

                virtual void setdPlasticMicroVelocityGradientdStateVariables( );

                virtual void setPreviousdPlasticMicroVelocityGradientdStateVariables( );

                virtual void setdPlasticGradientMicroVelocityGradientdMicroStress( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdMicroStress( );

                virtual void setdPlasticGradientMicroVelocityGradientdHigherOrderStress( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdHigherOrderStress( );

                virtual void setdPlasticGradientMicroVelocityGradientdF( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdF( );

                virtual void setdPlasticGradientMicroVelocityGradientdFn( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdFn( );

                virtual void setdPlasticGradientMicroVelocityGradientdChi( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdChi( );

                virtual void setdPlasticGradientMicroVelocityGradientdChin( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdChin( );

                virtual void setdPlasticGradientMicroVelocityGradientdGradChi( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdGradChi( );

                virtual void setdPlasticGradientMicroVelocityGradientdGradChin( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdGradChin( );

                virtual void setdPlasticGradientMicroVelocityGradientdStateVariables( );

                virtual void setPreviousdPlasticGradientMicroVelocityGradientdStateVariables( );

                virtual void setPlasticVelocityGradientsJacobians( const bool isPrevious );

                virtual void setPrecedingGradientMicroDeformation( );

                virtual void setPreviousPrecedingGradientMicroDeformation( );

                virtual void setPrecedingGradientMicroDeformation( const bool isPrevious );

                virtual void setdPrecedingGradientMicroDeformationdFn( );

                virtual void setdPrecedingGradientMicroDeformationdChi( );

                virtual void setdPrecedingGradientMicroDeformationdChin( );

                virtual void setdPrecedingGradientMicroDeformationdGradChi( );

                virtual void setdPrecedingGradientMicroDeformationdGradChin( );

                virtual void setPreviousdPrecedingGradientMicroDeformationdFn( );

                virtual void setPreviousdPrecedingGradientMicroDeformationdChi( );

                virtual void setPreviousdPrecedingGradientMicroDeformationdChin( );

                virtual void setPreviousdPrecedingGradientMicroDeformationdGradChi( );

                virtual void setPreviousdPrecedingGradientMicroDeformationdGradChin( );

                virtual void setPrecedingGradientMicroDeformationJacobians( const bool isPrevious );

                virtual void setUpdatedPlasticDeformationGradient( );

                virtual void setUpdatedPlasticMicroDeformation( );

                virtual void setUpdatedPlasticGradientMicroDeformation( );

                virtual void setPlasticDeformation( );

                virtual void setdUpdatedPlasticDeformationGradientdMacroStress( );

                virtual void setdUpdatedPlasticDeformationGradientdPreviousMacroStress( );

                virtual void setdUpdatedPlasticDeformationGradientdMicroStress( );

                virtual void setdUpdatedPlasticDeformationGradientdPreviousMicroStress( );

                virtual void setdUpdatedPlasticDeformationGradientdF( );

                virtual void setdUpdatedPlasticDeformationGradientdPreviousF( );

                virtual void setdUpdatedPlasticDeformationGradientdFn( );

                virtual void setdUpdatedPlasticDeformationGradientdPreviousFn( );

                virtual void setdUpdatedPlasticDeformationGradientdStateVariables( );

                virtual void setdUpdatedPlasticDeformationGradientdPreviousStateVariables( );

                virtual void setdUpdatedPlasticMicroDeformationdMicroStress( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousMicroStress( );

                virtual void setdUpdatedPlasticMicroDeformationdF( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousF( );

                virtual void setdUpdatedPlasticMicroDeformationdFn( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousFn( );

                virtual void setdUpdatedPlasticMicroDeformationdChi( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousChi( );

                virtual void setdUpdatedPlasticMicroDeformationdChin( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousChin( );

                virtual void setdUpdatedPlasticMicroDeformationdStateVariables( );

                virtual void setdUpdatedPlasticMicroDeformationdPreviousStateVariables( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdMacroStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousMacroStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdMicroStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousMicroStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdHigherOrderStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdF( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousF( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdFn( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousFn( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdChi( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousChi( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdChin( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousChin( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdGradChi( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousGradChi( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdGradChin( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousGradChin( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdStateVariables( );

                virtual void setdUpdatedPlasticGradientMicroDeformationdPreviousStateVariables( );

                virtual void setPlasticDeformationJacobians( const bool addPreviousGradients );

                virtual void setStateVariableResiduals( );

                virtual void setStateVariableJacobians( );

                virtual void setdStateVariableResidualsdD( );

                virtual void setdStateVariableResidualsdPreviousISVs( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdD( ) override;

            private:

                unsigned int _plasticConfigurationIndex; //! The index of the plastic configuration

                std::vector< unsigned int > _stateVariableIndices; //! The indices of the state variables in the global solve

                floatType _integrationParameter; //! The integration parameter (0 is explicit, 1 is implicit)

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroHardeningParameters,                             floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microHardeningParameters,                             floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientHardeningParameters,                     floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroFlowParameters,                                  floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microFlowParameters,                                  floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientFlowParameters,                          floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroYieldParameters,                                 floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microYieldParameters,                                 floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientYieldParameters,                         floatVector, unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroDrivingStress,                                   floatVector, setMacroDrivingStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, symmetricMicroDrivingStress,                          floatVector, setSymmetricMicroDrivingStress                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, higherOrderDrivingStress,                             floatVector, setHigherOrderDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdMacroStress,                      floatMatrix, setdMacroDrivingStressdMacroStress                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdMicroStress,             floatMatrix, setdSymmetricMicroDrivingStressdMicroStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdHigherOrderStress,          floatMatrix, setdHigherOrderDrivingStressdHigherOrderStress          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdF,                                floatMatrix, setdMacroDrivingStressdF                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdF,                       floatMatrix, setdSymmetricMicroDrivingStressdF                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdF,                          floatMatrix, setdHigherOrderDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChi,                        floatMatrix, setdHigherOrderDrivingStressdChi                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdFn,                               floatMatrix, setdMacroDrivingStressdFn                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdFn,                      floatMatrix, setdSymmetricMicroDrivingStressdFn                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdFn,                         floatMatrix, setdHigherOrderDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChin,                       floatMatrix, setdHigherOrderDrivingStressdChin                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroDrivingStress,                           floatVector, setPreviousMacroDrivingStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSymmetricMicroDrivingStress,                  floatVector, setPreviousSymmetricMicroDrivingStress                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHigherOrderDrivingStress,                     floatVector, setPreviousHigherOrderDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdMacroStress,              floatMatrix, setPreviousdMacroDrivingStressdMacroStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdMicroStress,     floatMatrix, setPreviousdSymmetricMicroDrivingStressdMicroStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdHigherOrderStress,  floatMatrix, setPreviousdHigherOrderDrivingStressdHigherOrderStress  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdF,                        floatMatrix, setPreviousdMacroDrivingStressdF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdF,               floatMatrix, setPreviousdSymmetricMicroDrivingStressdF               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdF,                  floatMatrix, setPreviousdHigherOrderDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdChi,                floatMatrix, setPreviousdHigherOrderDrivingStressdChi                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdFn,                       floatMatrix, setPreviousdMacroDrivingStressdFn                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdFn,              floatMatrix, setPreviousdSymmetricMicroDrivingStressdFn              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdFn,                 floatMatrix, setPreviousdHigherOrderDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdChin,               floatMatrix, setPreviousdHigherOrderDrivingStressdChin               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMultipliers,                                   floatVector, setPlasticMultipliers                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMultipliers,                           floatVector, setPreviousPlasticMultipliers                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStrainLikeISVs,                                floatVector, setPlasticStrainLikeISVs                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticStrainLikeISVs,                        floatVector, setPreviousPlasticStrainLikeISVs                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStateVariables,                                floatVector, setPlasticStateVariables                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticStateVariables,                        floatVector, setPreviousPlasticStateVariables                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroCohesion,                                        floatType,   setMacroCohesion                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroCohesiondStateVariables,                        floatVector, setdMacroCohesiondStateVariables                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microCohesion,                                        floatType,   setMicroCohesion                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroCohesiondStateVariables,                        floatVector, setdMicroCohesiondStateVariables                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microGradientCohesion,                                floatVector, setMicroGradientCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientCohesiondStateVariables,                floatMatrix, setdMicroGradientCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroCohesion,                                floatType,   setPreviousMacroCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroCohesiondStateVariables,                floatVector, setPreviousdMacroCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroCohesion,                                floatType,   setPreviousMicroCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroCohesiondStateVariables,                floatVector, setPreviousdMicroCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroGradientCohesion,                        floatVector, setPreviousMicroGradientCohesion                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientCohesiondStateVariables,        floatMatrix, setPreviousdMicroGradientCohesiondStateVariables        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroFlowdc,                                         floatType,   setdMacroFlowdc                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroFlowdDrivingStress,                             floatVector, setdMacroFlowdDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdStress,                     floatMatrix, setd2MacroFlowdDrivingStressdStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdF,                          floatMatrix, setd2MacroFlowdDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdFn,                         floatMatrix, setd2MacroFlowdDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroFlowdc,                                         floatType,   setdMicroFlowdc                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroFlowdDrivingStress,                             floatVector, setdMicroFlowdDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdStress,                     floatMatrix, setd2MicroFlowdDrivingStressdStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdF,                          floatMatrix, setd2MicroFlowdDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdFn,                         floatMatrix, setd2MicroFlowdDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientFlowdc,                                 floatMatrix, setdMicroGradientFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientFlowdDrivingStress,                     floatMatrix, setdMicroGradientFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdStress,             floatMatrix, setd2MicroGradientFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdF,                  floatMatrix, setd2MicroGradientFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdFn,                 floatMatrix, setd2MicroGradientFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdChi,                floatMatrix, setd2MicroGradientFlowdDrivingStressdChi                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdChin,               floatMatrix, setd2MicroGradientFlowdDrivingStressdChin               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroFlowdc,                                 floatType,   setPreviousdMacroFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroFlowdDrivingStress,                     floatVector, setPreviousdMacroFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdStress,             floatMatrix, setPreviousd2MacroFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdF,                  floatMatrix, setPreviousd2MacroFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdFn,                 floatMatrix, setPreviousd2MacroFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroFlowdc,                                 floatType,   setPreviousdMicroFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroFlowdDrivingStress,                     floatVector, setPreviousdMicroFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdStress,             floatMatrix, setPreviousd2MicroFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdF,                  floatMatrix, setPreviousd2MicroFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdFn,                 floatMatrix, setPreviousd2MicroFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientFlowdc,                         floatMatrix, setPreviousdMicroGradientFlowdc                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientFlowdDrivingStress,             floatMatrix, setPreviousdMicroGradientFlowdDrivingStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdStress,     floatMatrix, setPreviousd2MicroGradientFlowdDrivingStressdStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdF,          floatMatrix, setPreviousd2MicroGradientFlowdDrivingStressdF          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdFn,         floatMatrix, setPreviousd2MicroGradientFlowdDrivingStressdFn         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdChi,        floatMatrix, setPreviousd2MicroGradientFlowdDrivingStressdChi        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdChin,       floatMatrix, setPreviousd2MicroGradientFlowdDrivingStressdChin       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStrainLikeISVEvolutionRates,                   floatVector, setPlasticStrainLikeISVEvolutionRates                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStrainLikeISVEvolutionRatesdStateVariables,   floatMatrix, setdPlasticStrainLikeISVEvolutionRatesdStateVariables   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticStrainLikeISVEvolutionRates,           floatVector, setPreviousPlasticStrainLikeISVEvolutionRates           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticStrainLikeISVEvolutionRatesdStateVariables,  floatMatrix, setPreviousdPlasticStrainLikeISVEvolutionRatesdStateVariables  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticStrainLikeISVs,                         floatVector, setUpdatedPlasticStrainLikeISVs                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticStrainLikeISVsdStateVariables,         floatMatrix, setdUpdatedPlasticStrainLikeISVsdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticStrainLikeISVsdPreviousStateVariables, floatMatrix, setdUpdatedPlasticStrainLikeISVsdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroYield,                                           floatType,   setMacroYield                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddStress,                                   floatVector, setdMacroYielddStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddStateVariables,                           floatVector, setdMacroYielddStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddF,                                        floatVector, setdMacroYielddF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddFn,                                       floatVector, setdMacroYielddFn                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microYield,                                           floatType,   setMicroYield                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddStress,                                   floatVector, setdMicroYielddStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddStateVariables,                           floatVector, setdMicroYielddStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddF,                                        floatVector, setdMicroYielddF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddFn,                                       floatVector, setdMicroYielddFn                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microGradientYield,                                   floatVector, setMicroGradientYield                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddStress,                           floatMatrix, setdMicroGradientYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddStateVariables,                   floatMatrix, setdMicroGradientYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddF,                                floatMatrix, setdMicroGradientYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddFn,                               floatMatrix, setdMicroGradientYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddChi,                              floatMatrix, setdMicroGradientYielddChi                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddChin,                             floatMatrix, setdMicroGradientYielddChin                             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroYield,                                   floatType,   setPreviousMacroYield                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddStress,                           floatVector, setPreviousdMacroYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddStateVariables,                   floatVector, setPreviousdMacroYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddF,                                floatVector, setPreviousdMacroYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddFn,                               floatVector, setPreviousdMacroYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroYield,                                   floatType,   setPreviousMicroYield                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddStress,                           floatVector, setPreviousdMicroYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddStateVariables,                   floatVector, setPreviousdMicroYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddF,                                floatVector, setPreviousdMicroYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddFn,                               floatVector, setPreviousdMicroYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroGradientYield,                           floatVector, setPreviousMicroGradientYield                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddStress,                   floatMatrix, setPreviousdMicroGradientYielddStress                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddStateVariables,           floatMatrix, setPreviousdMicroGradientYielddStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddF,                        floatMatrix, setPreviousdMicroGradientYielddF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddFn,                       floatMatrix, setPreviousdMicroGradientYielddFn                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddChi,                      floatMatrix, setPreviousdMicroGradientYielddChi                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddChin,                     floatMatrix, setPreviousdMicroGradientYielddChin                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingDeformationGradient,                         floatVector, setPrecedingDeformationGradient                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdF,                      floatMatrix, setdPrecedingDeformationGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdFn,                     floatMatrix, setdPrecedingDeformationGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingDeformationGradient,                 floatVector, setPreviousPrecedingDeformationGradient                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingDeformationGradientdF,              floatMatrix, setPreviousdPrecedingDeformationGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingDeformationGradientdFn,             floatMatrix, setPreviousdPrecedingDeformationGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingMicroDeformation,                            floatVector, setPrecedingMicroDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingMicroDeformationdChi,                       floatMatrix, setdPrecedingMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingMicroDeformationdChin,                      floatMatrix, setdPrecedingMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingMicroDeformation,                    floatVector, setPreviousPrecedingMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingMicroDeformationdChi,               floatMatrix, setPreviousdPrecedingMicroDeformationdChi               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingMicroDeformationdChin,              floatMatrix, setPreviousdPrecedingMicroDeformationdChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMacroVelocityGradient,                         floatVector, setPlasticMacroVelocityGradient                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMacroVelocityGradient,                 floatVector, setPreviousPlasticMacroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdMacroStress,            floatMatrix, setdPlasticMacroVelocityGradientdMacroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdMacroStress,    floatMatrix, setPreviousdPlasticMacroVelocityGradientdMacroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdMicroStress,            floatMatrix, setdPlasticMacroVelocityGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdMicroStress,    floatMatrix, setPreviousdPlasticMacroVelocityGradientdMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdF,                      floatMatrix, setdPlasticMacroVelocityGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdF,              floatMatrix, setPreviousdPlasticMacroVelocityGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdFn,                     floatMatrix, setdPlasticMacroVelocityGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdFn,             floatMatrix, setPreviousdPlasticMacroVelocityGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdStateVariables,         floatMatrix, setdPlasticMacroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdStateVariables, floatMatrix, setPreviousdPlasticMacroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMicroVelocityGradient,                         floatVector, setPlasticMicroVelocityGradient                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMicroVelocityGradient,                 floatVector, setPreviousPlasticMicroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdMicroStress,            floatMatrix, setdPlasticMicroVelocityGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdMicroStress,    floatMatrix, setPreviousdPlasticMicroVelocityGradientdMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdF,                      floatMatrix, setdPlasticMicroVelocityGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdF,              floatMatrix, setPreviousdPlasticMicroVelocityGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdFn,                     floatMatrix, setdPlasticMicroVelocityGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdFn,             floatMatrix, setPreviousdPlasticMicroVelocityGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdChi,                    floatMatrix, setdPlasticMicroVelocityGradientdChi                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdChi,            floatMatrix, setPreviousdPlasticMicroVelocityGradientdChi            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdChin,                   floatMatrix, setdPlasticMicroVelocityGradientdChin                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdChin,           floatMatrix, setPreviousdPlasticMicroVelocityGradientdChin           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdStateVariables,         floatMatrix, setdPlasticMicroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdStateVariables, floatMatrix, setPreviousdPlasticMicroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingGradientMicroDeformation,                    floatVector, setPrecedingGradientMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdFn,                floatMatrix, setdPrecedingGradientMicroDeformationdFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdChi,               floatMatrix, setdPrecedingGradientMicroDeformationdChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdChin,              floatMatrix, setdPrecedingGradientMicroDeformationdChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdGradChi,           floatMatrix, setdPrecedingGradientMicroDeformationdGradChi           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdGradChin,          floatMatrix, setdPrecedingGradientMicroDeformationdGradChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingGradientMicroDeformation,            floatVector, setPreviousPrecedingGradientMicroDeformation            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdFn,        floatMatrix, setPreviousdPrecedingGradientMicroDeformationdFn        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdChi,       floatMatrix, setPreviousdPrecedingGradientMicroDeformationdChi       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdChin,      floatMatrix, setPreviousdPrecedingGradientMicroDeformationdChin      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdGradChi,   floatMatrix, setPreviousdPrecedingGradientMicroDeformationdGradChi   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdGradChin,  floatMatrix, setPreviousdPrecedingGradientMicroDeformationdGradChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticGradientMicroVelocityGradient,                 floatVector, setPlasticGradientMicroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticGradientMicroVelocityGradient,         floatVector, setPreviousPlasticGradientMicroVelocityGradient         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdMicroStress,         floatMatrix, setdPlasticGradientMicroVelocityGradientdMicroStress         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdMicroStress, floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdMicroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdHigherOrderStress,         floatMatrix, setdPlasticGradientMicroVelocityGradientdHigherOrderStress         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdHigherOrderStress, floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdHigherOrderStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdF,              floatMatrix, setdPlasticGradientMicroVelocityGradientdF             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdF,      floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdF     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdFn,             floatMatrix, setdPlasticGradientMicroVelocityGradientdFn            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdFn,     floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdFn    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdChi,            floatMatrix, setdPlasticGradientMicroVelocityGradientdChi           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdChi,    floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdChi   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdChin,           floatMatrix, setdPlasticGradientMicroVelocityGradientdChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdChin,   floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdGradChi,            floatMatrix, setdPlasticGradientMicroVelocityGradientdGradChi           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdGradChi,    floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdGradChi   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdGradChin,           floatMatrix, setdPlasticGradientMicroVelocityGradientdGradChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdGradChin,   floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdGradChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdStateVariables,           floatMatrix, setdPlasticGradientMicroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdStateVariables,   floatMatrix, setPreviousdPlasticGradientMicroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticDeformationGradient,                         floatVector, setUpdatedPlasticDeformationGradient                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticMicroDeformation,                            floatVector, setUpdatedPlasticMicroDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticGradientMicroDeformation,                    floatVector, setUpdatedPlasticGradientMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdMacroStress,            floatMatrix, setdUpdatedPlasticDeformationGradientdMacroStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousMacroStress,    floatMatrix, setdUpdatedPlasticDeformationGradientdPreviousMacroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdMicroStress,            floatMatrix, setdUpdatedPlasticDeformationGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousMicroStress,    floatMatrix, setdUpdatedPlasticDeformationGradientdPreviousMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdF,                      floatMatrix, setdUpdatedPlasticDeformationGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousF,              floatMatrix, setdUpdatedPlasticDeformationGradientdPreviousF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdFn,                     floatMatrix, setdUpdatedPlasticDeformationGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousFn,             floatMatrix, setdUpdatedPlasticDeformationGradientdPreviousFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdStateVariables,         floatMatrix, setdUpdatedPlasticDeformationGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousStateVariables, floatMatrix, setdUpdatedPlasticDeformationGradientdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdMicroStress,               floatMatrix, setdUpdatedPlasticMicroDeformationdMicroStress               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousMicroStress,       floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousMicroStress       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdF,                         floatMatrix, setdUpdatedPlasticMicroDeformationdF                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousF,                 floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdFn,                        floatMatrix, setdUpdatedPlasticMicroDeformationdFn                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousFn,                floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdChi,                       floatMatrix, setdUpdatedPlasticMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousChi,               floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdChin,                      floatMatrix, setdUpdatedPlasticMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousChin,              floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdStateVariables,            floatMatrix, setdUpdatedPlasticMicroDeformationdStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousStateVariables,    floatMatrix, setdUpdatedPlasticMicroDeformationdPreviousStateVariables    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdMacroStress,         floatMatrix, setdUpdatedPlasticGradientMicroDeformationdMacroStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress, floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousMacroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdMicroStress,         floatMatrix, setdUpdatedPlasticGradientMicroDeformationdMicroStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress, floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousMicroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdHigherOrderStress,         floatMatrix, setdUpdatedPlasticGradientMicroDeformationdHigherOrderStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress, floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdF,                         floatMatrix, setdUpdatedPlasticGradientMicroDeformationdF                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousF,                 floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdFn,                        floatMatrix, setdUpdatedPlasticGradientMicroDeformationdFn                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousFn,                floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdChi,                       floatMatrix, setdUpdatedPlasticGradientMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousChi,               floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdChin,                      floatMatrix, setdUpdatedPlasticGradientMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousChin,              floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdGradChi,                   floatMatrix, setdUpdatedPlasticGradientMicroDeformationdGradChi                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousGradChi,           floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousGradChi           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdGradChin,                  floatMatrix, setdUpdatedPlasticGradientMicroDeformationdGradChin                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousGradChin,          floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousGradChin          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdStateVariables,            floatMatrix, setdUpdatedPlasticGradientMicroDeformationdStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables,    floatMatrix, setdUpdatedPlasticGradientMicroDeformationdPreviousStateVariables    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableResiduals,                                            floatVector, setStateVariableResiduals                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableJacobians,                                            floatMatrix, setStateVariableJacobians                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableResidualsdD,                                         floatMatrix, setdStateVariableResidualsdD                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableResidualsdPreviousISVs,                              floatMatrix, setdStateVariableResidualsdPreviousISVs                              )

        };

    }

}

#endif
