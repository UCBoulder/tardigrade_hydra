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

        typedef std::vector< floatType > seventhOrderTensor; //!< A seventh order tensor
        typedef std::vector< floatType > eighthOrderTensor; //!< A eighth order tensor

        template <typename T> int sgn(T val){
            return (T(0) < val) - (val < T(0));
        }

        void computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B );

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &stressMeasure,
                                                           const floatType         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue );

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &stressMeasure,
                                                           const floatType         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue, secondOrderTensor &dFdStress, floatType &dFdc,
                                                           secondOrderTensor &dFdPrecedingF, double tol = 1e-9 );

        void computeSecondOrderDruckerPragerYieldEquation( const secondOrderTensor &stressMeasure,
                                                           const floatType         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           floatType &yieldValue, secondOrderTensor &dFdStress, floatType &dFdc,
                                                           secondOrderTensor &dFdPrecedingF, fourthOrderTensor &d2FdStress2,
                                                           fourthOrderTensor &d2FdStressdPrecedingF, double tol = 1e-9 );

        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor  &stressMeasure,
                                                           const dimVector         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           dimVector &yieldValue );
    
        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor  &stressMeasure,
                                                           const dimVector         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           dimVector &yieldValue, fourthOrderTensor &dFdStress, dimVector &dFdc,
                                                           thirdOrderTensor &dFdPrecedingF );
    
        void computeHigherOrderDruckerPragerYieldEquation( const thirdOrderTensor  &stressMeasure,
                                                           const dimVector         &cohesion,
                                                           const secondOrderTensor &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           dimVector &yieldValue, thirdOrderTensor &dFdStress, dimVector &dFdc,
                                                           thirdOrderTensor &dFdPrecedingF, seventhOrderTensor &d2FdStress2,
                                                           sixthOrderTensor &d2FdStressdPrecedingF );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient,
                                                  secondOrderTensor &dPlasticMacroLdMacroGamma,
                                                  secondOrderTensor &dPlasticMacroLdMicroGamma );

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const secondOrderTensor &inverseElasticRightCauchyGreen,
                                                  const secondOrderTensor &macroFlowDirection,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMacroVelocityGradient,
                                                  secondOrderTensor &dPlasticMacroLdMacroGamma,
                                                  secondOrderTensor &dPlasticMacroLdMicroGamma,
                                                  fourthOrderTensor &dPlasticMacroLdElasticRCG,
                                                  fourthOrderTensor &dPlasticMacroLdMacroFlowDirection,
                                                  fourthOrderTensor &dPlasticMacroLdMicroFlowDirection );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient,
                                                  secondOrderTensor &dPlasticMicroLdMicroGamma );

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const secondOrderTensor &elasticMicroRightCauchyGreen,
                                                  const secondOrderTensor &elasticPsi, const secondOrderTensor &inverseElasticPsi,
                                                  const secondOrderTensor &microFlowDirection,
                                                  secondOrderTensor &plasticMicroVelocityGradient,
                                                  secondOrderTensor &dPlasticMicroLdMicroGamma,
                                                  fourthOrderTensor &dPlasticMicroLdElasticMicroRCG,
                                                  fourthOrderTensor &dPlasticMicroLdElasticPsi,
                                                  fourthOrderTensor &dPlasticMicroLdMicroFlowDirection );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          thirdOrderTensor &skewTerm );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          thirdOrderTensor &skewTerm,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL,
                                                          fifthOrderTensor &dPlasticMicroGradientLdElasticPsi,
                                                          sixthOrderTensor &dPlasticMicroGradientLdElasticGamma,
                                                          sixthOrderTensor &dPlasticMicroGradientLdMicroGradientFlowDirection );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          thirdOrderTensor &skewTerm,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL );

        void computePlasticMicroGradientVelocityGradient( const dimVector &microGradientGamma, const secondOrderTensor &elasticPsi,
                                                          const secondOrderTensor &inverseElasticPsi, const thirdOrderTensor &elasticGamma,
                                                          const thirdOrderTensor &microGradientFlowDirection,
                                                          const secondOrderTensor &plasticMicroVelocityGradient,
                                                          thirdOrderTensor &plasticMicroGradientVelocityGradient,
                                                          fourthOrderTensor &dPlasticMicroGradientLdMicroGradientGamma,
                                                          fifthOrderTensor &dPlasticMicroGradientLdPlasticMicroL,
                                                          fifthOrderTensor &dPlasticMicroGradientLdElasticPsi,
                                                          sixthOrderTensor &dPlasticMicroGradientLdElasticGamma,
                                                          sixthOrderTensor &dPlasticMicroGradientLdMicroGradientFlowDirection );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        sixthOrderTensor &LHS,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticMicroGradChi( const variableType &Dt,
                                        const secondOrderTensor &currentPlasticMicroDeformation,
                                        const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                        const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroDeformation,
                                        const thirdOrderTensor &previousPlasticMicroGradient,
                                        const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                        const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                        const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                        thirdOrderTensor &currentPlasticMicroGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                                        fifthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                                        sixthOrderTensor &dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor &currentPlasticMicroGradient,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor &currentPlasticMicroGradient,
                                       fourthOrderTensor &dPlasticFdPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPlasticMicroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPlasticMacroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPlasticMicroL,
                                       sixthOrderTensor &dPlasticMicroGradientdPlasticMicroGradientL,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        void evolvePlasticDeformation( const variableType &Dt,
                                       const secondOrderTensor &currentPlasticMacroVelocityGradient,
                                       const secondOrderTensor &currentPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &currentPlasticMicroGradientVelocityGradient,
                                       const secondOrderTensor &previousPlasticDeformationGradient,
                                       const secondOrderTensor &previousPlasticMicroDeformation,
                                       const thirdOrderTensor &previousPlasticMicroGradient,
                                       const secondOrderTensor &previousPlasticMacroVelocityGradient,
                                       const secondOrderTensor &previousPlasticMicroVelocityGradient,
                                       const thirdOrderTensor &previousPlasticMicroGradientVelocityGradient,
                                       secondOrderTensor &currentPlasticDeformationGradient,
                                       secondOrderTensor &currentPlasticMicroDeformation,
                                       thirdOrderTensor &currentPlasticMicroGradient,
                                       fourthOrderTensor &dPlasticFdPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPlasticMicroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPlasticMacroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPlasticMicroL,
                                       sixthOrderTensor &dPlasticMicroGradientdPlasticMicroGradientL,
                                       fourthOrderTensor &dPlasticFdPreviousPlasticF,
                                       fourthOrderTensor &dPlasticFdPreviousPlasticMacroL,
                                       fourthOrderTensor &dPlasticMicroDeformationdPreviousPlasticMicroDeformation,
                                       fourthOrderTensor &dPlasticMicroDeformationdPreviousPlasticMicroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                       sixthOrderTensor &dPlasticMicroGradientdPreviousPlasticMicroGradient,
                                       fifthOrderTensor &dPlasticMicroGradientdPreviousPlasticMacroL,
                                       fifthOrderTensor &dPlasticMicroGradientdPreviousPlasticMicroL,
                                       sixthOrderTensor &dPlasticMicroGradientdPreviousPlasticMicroGradientL,
                                       const parameterType alphaMacro = 0.5,
                                       const parameterType alphaMicro = 0.5,
                                       const parameterType alphaMicroGradient = 0.5 );

        variableType weakMac( const variableType &x, const variableType &a );

        variableType weakMac( const variableType &x, const variableType &a, variableType &dmacdx );

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
                          const std::vector< unsigned int > &stateVariableIndices, const floatVector &parameters, const floatType integrationParameter = 0.5,
                          const bool useWeakenedMacaulay = false, const floatType weakenedMacaulayParameter=10, const floatType plasticMultiplierBarrierModulus=1000. )
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
                     * \param &useWeakenedMacaulay: A flag for whether to use a weakened Macaulay bracket or not (defaults to false)
                     * \param &weakenedMacaulayParameter: The value of the parameter for the weakened Macaulay bracket (defaults to 10)
                     * \param &plasticMultiplierBarrierModulus: The barrier modulus to make sure that the plastic multipliers never go negative (defaults to 1000)
                     */

                    _plasticConfigurationIndex = plasticConfigurationIndex;

                    _numPlasticMultipliers = 5; //Hard-coding this because we couldn't handle a different number anyway.

                    _numPlasticStrainLikeStateVariables = 5; //Hard-coding this because we couldn't handle a different number anyway.

                    _stateVariableIndices = stateVariableIndices;

                    _integrationParameter = integrationParameter;

                    _useWeakenedMacaulay = useWeakenedMacaulay;

                    _weakenedMacaulayParameter = weakenedMacaulayParameter;

                    _plasticMultiplierBarrierModulus = plasticMultiplierBarrierModulus;

                    TARDIGRADE_ERROR_TOOLS_CATCH( extractMaterialParameters( parameters ) );

                }

                //!Get the number of plastic multipliers expected for the problem
                const unsigned int* getNumPlasticMultipliers( ){ return &_numPlasticMultipliers; }

                //!Get the number of plastic multipliers expected for the problem
                const unsigned int* getNumStrainLikePlasticStateVariables( ){ return &_numPlasticStrainLikeStateVariables; }

                const unsigned int* getPlasticConfigurationIndex( );

                const std::vector< unsigned int >* getStateVariableIndices( );

                const floatType* getIntegrationParameter( );

                //!Return the flag for whether to use the weakened Macaulay bracket or not
                const bool *useWeakenedMacaulay( ){ return &_useWeakenedMacaulay; }

                //!Return the value of the Weakened Macaulay parameter
                const floatType *getWeakenedMacaulayParameter( ){ return &_weakenedMacaulayParameter; }

                //!Return the value of the barrier modulus to prevent the plastic multipliers from becoming negative
                const floatType *getPlasticMultiplierBarrierModulus( ){ return &_plasticMultiplierBarrierModulus; }

                //!Return the value of the consistency condition modulus
                const floatType *getConsistencyConditionModulus( ){ return &_consistencyConditionModulus; }

                virtual void projectSuggestedX( std::vector< floatType > &trialX,
                                                const std::vector< floatType > &Xp ) override;

                //!Get the maximum allowable value for the norm of the change in macro plastic deformation for a single nonlinear increment
                const floatType* getMaxMacroPlasticDeltaNorm( ){ return &_maxMacroPlasticDeltaNorm; }

                //!Get the maximum allowable value for the norm of the change in micro plastic deformation for a single nonlinear increment
                const floatType* getMaxMicroPlasticDeltaNorm( ){ return &_maxMicroPlasticDeltaNorm; }

                //!Get the maximum allowable value for the norm of the change in micro gradient plastic deformation for a single nonlinear increment
                const floatType* getMaxMicroGradientPlasticDeltaNorm( ){ return &_maxMicroGradientPlasticDeltaNorm; }

                void setMaxMacroPlasticDeltaNorm( const floatType &value ){
                    /*!
                     * Set the maximum allowable value for the norm of the change in macro plastic deformation for a single nonlinear increment
                     *
                     * \param &value: The parameter value
                     */

                    _maxMacroPlasticDeltaNorm = value;

                }

                void setMaxMicroPlasticDeltaNorm( const floatType &value ){
                    /*!
                     * Set the maximum allowable value for the norm of the change in micro plastic deformation for a single nonlinear increment
                     *
                     * \param &value: The parameter value
                     */

                    _maxMicroPlasticDeltaNorm = value;

                }

                void setMaxMicroGradientPlasticDeltaNorm( const floatType &value ){
                    /*!
                     * Set the maximum allowable value for the norm of the change in micro gradient plastic deformation for a single nonlinear increment
                     *
                     * \param &value: The parameter value
                     */

                    _maxMicroGradientPlasticDeltaNorm = value;

                }

//                virtual void suggestInitialIterateValues( std::vector< unsigned int >   &indices,
//                                                          std::vector< floatType > &values ) override;

            protected:

                bool _useWeakenedMacaulay; //!< Flag for whether to use the weak Macaulay brackets or not

                floatType _weakenedMacaulayParameter; //!< The weakening parameter for the weak Macaulay brackets

                floatType _plasticMultiplierBarrierModulus; //!< The modulus applied to the plastic multipliers to make sure they are never negative

                floatType _consistencyConditionModulus = 0.1; //!< Modulus for the consistency condition

                unsigned int _numPlasticMultipliers; //!< The number of plastic multipliers. Hard coded to 5 but setting as a variable just in case

                unsigned int _numPlasticStrainLikeStateVariables = 5; //The number of plastic strain-like state variables. Hard coded to 5 but setting as a variable just in case

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

                virtual void setdRdT( ) override;

                const floatType *getMinCohesion( ){ return &_minCohesion; }

                const floatType *getSmoothRatio( ){ return &_smoothRatio; };

                virtual double softLinearCohesion( const floatType &Z, const floatType &A, const floatType &c0, const floatType &rc, const floatType &cf );

                virtual double softLinearCohesionDerivative( const floatType &Z, const floatType &A, const floatType &c0, const floatType &rc, const floatType &cf );

            private:

                unsigned int _plasticConfigurationIndex; //! The index of the plastic configuration

                std::vector< unsigned int > _stateVariableIndices; //! The indices of the state variables in the global solve

                floatType _integrationParameter; //! The integration parameter (0 is explicit, 1 is implicit)

                floatType _maxMacroPlasticDeltaNorm = 1.; //!< The maximum allowable value of the norm of the change in macro plasticity for a given nonlinear iteration

                floatType _maxMicroPlasticDeltaNorm = 1.; //!< The maximum allowable value of the norm of the change in micro plasticity for a given nonlinear iteration

                floatType _maxMicroGradientPlasticDeltaNorm = 1.; //!< The maximum allowable value of the norm of the change in micro gradient plasticity for a given nonlinear iteration

                floatType _smoothRatio = 0.1; //!< The ratio of the initial cohesion value at which to start smoothing the response to zero

                floatType _minCohesion = 1e-2; //!< The minimum allowable value of the cohesion

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroHardeningParameters,                             floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microHardeningParameters,                             floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientHardeningParameters,                     floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroFlowParameters,                                  floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microFlowParameters,                                  floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientFlowParameters,                          floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, macroYieldParameters,                                 floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microYieldParameters,                                 floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(  private, microGradientYieldParameters,                         floatVector,       unexpectedError                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroDrivingStress,                                   secondOrderTensor, setMacroDrivingStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, symmetricMicroDrivingStress,                          secondOrderTensor, setSymmetricMicroDrivingStress                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, higherOrderDrivingStress,                             thirdOrderTensor,  setHigherOrderDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdMacroStress,                      fourthOrderTensor, setdMacroDrivingStressdMacroStress                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdMicroStress,             fourthOrderTensor, setdSymmetricMicroDrivingStressdMicroStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdHigherOrderStress,          sixthOrderTensor,  setdHigherOrderDrivingStressdHigherOrderStress          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdF,                                fourthOrderTensor, setdMacroDrivingStressdF                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdF,                       fourthOrderTensor, setdSymmetricMicroDrivingStressdF                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdF,                          fifthOrderTensor,  setdHigherOrderDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChi,                        fifthOrderTensor,  setdHigherOrderDrivingStressdChi                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroDrivingStressdFn,                               floatVector, setdMacroDrivingStressdFn                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dSymmetricMicroDrivingStressdFn,                      floatVector, setdSymmetricMicroDrivingStressdFn                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdFn,                         floatVector, setdHigherOrderDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dHigherOrderDrivingStressdChin,                       floatVector, setdHigherOrderDrivingStressdChin                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroDrivingStress,                           secondOrderTensor, setPreviousMacroDrivingStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousSymmetricMicroDrivingStress,                  secondOrderTensor, setPreviousSymmetricMicroDrivingStress                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousHigherOrderDrivingStress,                     thirdOrderTensor,  setPreviousHigherOrderDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdMacroStress,              fourthOrderTensor, setPreviousdMacroDrivingStressdMacroStress              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdMicroStress,     fourthOrderTensor, setPreviousdSymmetricMicroDrivingStressdMicroStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdHigherOrderStress,  sixthOrderTensor,  setPreviousdHigherOrderDrivingStressdHigherOrderStress  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdF,                        fourthOrderTensor, setPreviousdMacroDrivingStressdF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdF,               fourthOrderTensor, setPreviousdSymmetricMicroDrivingStressdF               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdF,                  fifthOrderTensor,  setPreviousdHigherOrderDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdChi,                fifthOrderTensor,  setPreviousdHigherOrderDrivingStressdChi                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroDrivingStressdFn,                       floatVector, setPreviousdMacroDrivingStressdFn                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdSymmetricMicroDrivingStressdFn,              floatVector, setPreviousdSymmetricMicroDrivingStressdFn              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdFn,                 floatVector, setPreviousdHigherOrderDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdHigherOrderDrivingStressdChin,               floatVector, setPreviousdHigherOrderDrivingStressdChin               )

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

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microGradientCohesion,                                dimVector,   setMicroGradientCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientCohesiondStateVariables,                floatVector, setdMicroGradientCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroCohesion,                                floatType,   setPreviousMacroCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroCohesiondStateVariables,                floatVector, setPreviousdMacroCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroCohesion,                                floatType,   setPreviousMicroCohesion                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroCohesiondStateVariables,                floatVector, setPreviousdMicroCohesiondStateVariables                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroGradientCohesion,                        dimVector,   setPreviousMicroGradientCohesion                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientCohesiondStateVariables,        floatVector, setPreviousdMicroGradientCohesiondStateVariables        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroFlowdc,                                         floatType,   setdMacroFlowdc                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroFlowdDrivingStress,                             secondOrderTensor, setdMacroFlowdDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdStress,                     fourthOrderTensor, setd2MacroFlowdDrivingStressdStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdF,                          fourthOrderTensor, setd2MacroFlowdDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MacroFlowdDrivingStressdFn,                         floatVector, setd2MacroFlowdDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroFlowdc,                                         floatType,   setdMicroFlowdc                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroFlowdDrivingStress,                             secondOrderTensor, setdMicroFlowdDrivingStress                             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdStress,                     fourthOrderTensor, setd2MicroFlowdDrivingStressdStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdF,                          fourthOrderTensor, setd2MicroFlowdDrivingStressdF                          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroFlowdDrivingStressdFn,                         floatVector, setd2MicroFlowdDrivingStressdFn                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientFlowdc,                                 secondOrderTensor, setdMicroGradientFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientFlowdDrivingStress,                     fourthOrderTensor, setdMicroGradientFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdStress,             fifthOrderTensor,  setd2MicroGradientFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdF,                  fifthOrderTensor,  setd2MicroGradientFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdFn,                 floatVector, setd2MicroGradientFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdChi,                fifthOrderTensor,  setd2MicroGradientFlowdDrivingStressdChi                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, d2MicroGradientFlowdDrivingStressdChin,               floatVector, setd2MicroGradientFlowdDrivingStressdChin               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroFlowdc,                                 floatType,         setPreviousdMacroFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroFlowdDrivingStress,                     secondOrderTensor, setPreviousdMacroFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdStress,             fourthOrderTensor, setPreviousd2MacroFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdF,                  fourthOrderTensor, setPreviousd2MacroFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MacroFlowdDrivingStressdFn,                 floatVector, setPreviousd2MacroFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroFlowdc,                                 floatType,   setPreviousdMicroFlowdc                                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroFlowdDrivingStress,                     secondOrderTensor, setPreviousdMicroFlowdDrivingStress                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdStress,             fourthOrderTensor, setPreviousd2MicroFlowdDrivingStressdStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdF,                  fourthOrderTensor, setPreviousd2MicroFlowdDrivingStressdF                  )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroFlowdDrivingStressdFn,                 floatVector, setPreviousd2MicroFlowdDrivingStressdFn                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientFlowdc,                         secondOrderTensor,  setPreviousdMicroGradientFlowdc                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientFlowdDrivingStress,             fourthOrderTensor,  setPreviousdMicroGradientFlowdDrivingStress             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdStress,     seventhOrderTensor, setPreviousd2MicroGradientFlowdDrivingStressdStress     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdF,          sixthOrderTensor,   setPreviousd2MicroGradientFlowdDrivingStressdF          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdFn,         floatVector, setPreviousd2MicroGradientFlowdDrivingStressdFn         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdChi,        sixthOrderTensor,   setPreviousd2MicroGradientFlowdDrivingStressdChi        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousd2MicroGradientFlowdDrivingStressdChin,       floatVector, setPreviousd2MicroGradientFlowdDrivingStressdChin       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticStrainLikeISVEvolutionRates,                   floatVector, setPlasticStrainLikeISVEvolutionRates                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticStrainLikeISVEvolutionRatesdStateVariables,   floatVector, setdPlasticStrainLikeISVEvolutionRatesdStateVariables   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticStrainLikeISVEvolutionRates,           floatVector, setPreviousPlasticStrainLikeISVEvolutionRates           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticStrainLikeISVEvolutionRatesdStateVariables,  floatVector, setPreviousdPlasticStrainLikeISVEvolutionRatesdStateVariables  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticStrainLikeISVs,                         floatVector, setUpdatedPlasticStrainLikeISVs                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticStrainLikeISVsdStateVariables,         floatVector, setdUpdatedPlasticStrainLikeISVsdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticStrainLikeISVsdPreviousStateVariables, floatVector, setdUpdatedPlasticStrainLikeISVsdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, macroYield,                                           floatType,   setMacroYield                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddStress,                                   secondOrderTensor, setdMacroYielddStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddStateVariables,                           floatVector, setdMacroYielddStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddF,                                        secondOrderTensor, setdMacroYielddF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMacroYielddFn,                                       floatVector, setdMacroYielddFn                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microYield,                                           floatType,   setMicroYield                                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddStress,                                   secondOrderTensor, setdMicroYielddStress                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddStateVariables,                           floatVector, setdMicroYielddStateVariables                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddF,                                        secondOrderTensor, setdMicroYielddF                                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroYielddFn,                                       floatVector, setdMicroYielddFn                                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microGradientYield,                                   floatVector, setMicroGradientYield                                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddStress,                           fourthOrderTensor, setdMicroGradientYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddStateVariables,                   floatVector, setdMicroGradientYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddF,                                fourthOrderTensor, setdMicroGradientYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddFn,                               floatVector, setdMicroGradientYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddChi,                              fourthOrderTensor, setdMicroGradientYielddChi                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dMicroGradientYielddChin,                             floatVector, setdMicroGradientYielddChin                             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMacroYield,                                   floatType,   setPreviousMacroYield                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddStress,                           secondOrderTensor, setPreviousdMacroYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddStateVariables,                   floatVector, setPreviousdMacroYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddF,                                secondOrderTensor, setPreviousdMacroYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMacroYielddFn,                               floatVector, setPreviousdMacroYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroYield,                                   floatType,   setPreviousMicroYield                                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddStress,                           secondOrderTensor, setPreviousdMicroYielddStress                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddStateVariables,                   floatVector, setPreviousdMicroYielddStateVariables                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddF,                                secondOrderTensor, setPreviousdMicroYielddF                                )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroYielddFn,                               floatVector, setPreviousdMicroYielddFn                               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroGradientYield,                           dimVector,   setPreviousMicroGradientYield                           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddStress,                   fourthOrderTensor, setPreviousdMicroGradientYielddStress                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddStateVariables,           floatVector, setPreviousdMicroGradientYielddStateVariables           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddF,                        fourthOrderTensor, setPreviousdMicroGradientYielddF                        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddFn,                       floatVector, setPreviousdMicroGradientYielddFn                       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddChi,                      fourthOrderTensor, setPreviousdMicroGradientYielddChi                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdMicroGradientYielddChin,                     floatVector, setPreviousdMicroGradientYielddChin                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingDeformationGradient,                         secondOrderTensor, setPrecedingDeformationGradient                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdF,                      fourthOrderTensor, setdPrecedingDeformationGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingDeformationGradientdFn,                     floatVector, setdPrecedingDeformationGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingDeformationGradient,                 secondOrderTensor, setPreviousPrecedingDeformationGradient                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingDeformationGradientdF,              fourthOrderTensor, setPreviousdPrecedingDeformationGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingDeformationGradientdFn,             floatVector, setPreviousdPrecedingDeformationGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingMicroDeformation,                            secondOrderTensor, setPrecedingMicroDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingMicroDeformationdChi,                       fourthOrderTensor, setdPrecedingMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingMicroDeformationdChin,                      floatVector, setdPrecedingMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingMicroDeformation,                    secondOrderTensor, setPreviousPrecedingMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingMicroDeformationdChi,               fourthOrderTensor, setPreviousdPrecedingMicroDeformationdChi               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingMicroDeformationdChin,              floatVector, setPreviousdPrecedingMicroDeformationdChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMacroVelocityGradient,                         secondOrderTensor, setPlasticMacroVelocityGradient                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMacroVelocityGradient,                 secondOrderTensor, setPreviousPlasticMacroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdMacroStress,            fourthOrderTensor, setdPlasticMacroVelocityGradientdMacroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdMacroStress,    fourthOrderTensor, setPreviousdPlasticMacroVelocityGradientdMacroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdMicroStress,            fourthOrderTensor, setdPlasticMacroVelocityGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdMicroStress,    fourthOrderTensor, setPreviousdPlasticMacroVelocityGradientdMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdF,                      fourthOrderTensor, setdPlasticMacroVelocityGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdF,              fourthOrderTensor, setPreviousdPlasticMacroVelocityGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdFn,                     floatVector, setdPlasticMacroVelocityGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdFn,             floatVector, setPreviousdPlasticMacroVelocityGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMacroVelocityGradientdStateVariables,         floatVector, setdPlasticMacroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMacroVelocityGradientdStateVariables, floatVector, setPreviousdPlasticMacroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticMicroVelocityGradient,                         secondOrderTensor, setPlasticMicroVelocityGradient                         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticMicroVelocityGradient,                 secondOrderTensor, setPreviousPlasticMicroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdMicroStress,            fourthOrderTensor, setdPlasticMicroVelocityGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdMicroStress,    fourthOrderTensor, setPreviousdPlasticMicroVelocityGradientdMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdF,                      fourthOrderTensor, setdPlasticMicroVelocityGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdF,              fourthOrderTensor, setPreviousdPlasticMicroVelocityGradientdF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdFn,                     floatVector, setdPlasticMicroVelocityGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdFn,             floatVector, setPreviousdPlasticMicroVelocityGradientdFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdChi,                    fourthOrderTensor, setdPlasticMicroVelocityGradientdChi                    )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdChi,            fourthOrderTensor, setPreviousdPlasticMicroVelocityGradientdChi            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdChin,                   floatVector, setdPlasticMicroVelocityGradientdChin                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdChin,           floatVector, setPreviousdPlasticMicroVelocityGradientdChin           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticMicroVelocityGradientdStateVariables,         floatVector, setdPlasticMicroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticMicroVelocityGradientdStateVariables, floatVector, setPreviousdPlasticMicroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, precedingGradientMicroDeformation,                    thirdOrderTensor, setPrecedingGradientMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdFn,                floatVector, setdPrecedingGradientMicroDeformationdFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdChi,               fifthOrderTensor, setdPrecedingGradientMicroDeformationdChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdChin,              floatVector, setdPrecedingGradientMicroDeformationdChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdGradChi,           sixthOrderTensor, setdPrecedingGradientMicroDeformationdGradChi           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPrecedingGradientMicroDeformationdGradChin,          floatVector, setdPrecedingGradientMicroDeformationdGradChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPrecedingGradientMicroDeformation,            thirdOrderTensor, setPreviousPrecedingGradientMicroDeformation            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdFn,        floatVector, setPreviousdPrecedingGradientMicroDeformationdFn        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdChi,       fifthOrderTensor, setPreviousdPrecedingGradientMicroDeformationdChi       )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdChin,      floatVector, setPreviousdPrecedingGradientMicroDeformationdChin      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdGradChi,   sixthOrderTensor, setPreviousdPrecedingGradientMicroDeformationdGradChi   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPrecedingGradientMicroDeformationdGradChin,  floatVector, setPreviousdPrecedingGradientMicroDeformationdGradChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, plasticGradientMicroVelocityGradient,                 thirdOrderTensor, setPlasticGradientMicroVelocityGradient                 )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPlasticGradientMicroVelocityGradient,         thirdOrderTensor, setPreviousPlasticGradientMicroVelocityGradient         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdMicroStress,         fifthOrderTensor, setdPlasticGradientMicroVelocityGradientdMicroStress         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdMicroStress, fifthOrderTensor, setPreviousdPlasticGradientMicroVelocityGradientdMicroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdHigherOrderStress,         sixthOrderTensor, setdPlasticGradientMicroVelocityGradientdHigherOrderStress         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdHigherOrderStress, sixthOrderTensor, setPreviousdPlasticGradientMicroVelocityGradientdHigherOrderStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdF,              fifthOrderTensor, setdPlasticGradientMicroVelocityGradientdF             )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdF,      fifthOrderTensor, setPreviousdPlasticGradientMicroVelocityGradientdF     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdFn,             floatVector, setdPlasticGradientMicroVelocityGradientdFn            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdFn,     floatVector, setPreviousdPlasticGradientMicroVelocityGradientdFn    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdChi,            fifthOrderTensor, setdPlasticGradientMicroVelocityGradientdChi           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdChi,    fifthOrderTensor, setPreviousdPlasticGradientMicroVelocityGradientdChi   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdChin,           floatVector, setdPlasticGradientMicroVelocityGradientdChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdChin,   floatVector, setPreviousdPlasticGradientMicroVelocityGradientdChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdGradChi,            sixthOrderTensor, setdPlasticGradientMicroVelocityGradientdGradChi           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdGradChi,    sixthOrderTensor, setPreviousdPlasticGradientMicroVelocityGradientdGradChi   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdGradChin,           floatVector, setdPlasticGradientMicroVelocityGradientdGradChin          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdGradChin,   floatVector, setPreviousdPlasticGradientMicroVelocityGradientdGradChin  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPlasticGradientMicroVelocityGradientdStateVariables,           floatVector, setdPlasticGradientMicroVelocityGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPlasticGradientMicroVelocityGradientdStateVariables,   floatVector, setPreviousdPlasticGradientMicroVelocityGradientdStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticDeformationGradient,                         secondOrderTensor, setUpdatedPlasticDeformationGradient                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticMicroDeformation,                            secondOrderTensor, setUpdatedPlasticMicroDeformation                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, updatedPlasticGradientMicroDeformation,                    thirdOrderTensor, setUpdatedPlasticGradientMicroDeformation                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdMacroStress,            fourthOrderTensor, setdUpdatedPlasticDeformationGradientdMacroStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousMacroStress,    fourthOrderTensor, setdUpdatedPlasticDeformationGradientdPreviousMacroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdMicroStress,            fourthOrderTensor, setdUpdatedPlasticDeformationGradientdMicroStress            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousMicroStress,    fourthOrderTensor, setdUpdatedPlasticDeformationGradientdPreviousMicroStress    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdF,                      fourthOrderTensor, setdUpdatedPlasticDeformationGradientdF                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousF,              fourthOrderTensor, setdUpdatedPlasticDeformationGradientdPreviousF              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdFn,                     floatVector, setdUpdatedPlasticDeformationGradientdFn                     )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousFn,             floatVector, setdUpdatedPlasticDeformationGradientdPreviousFn             )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdStateVariables,         floatVector, setdUpdatedPlasticDeformationGradientdStateVariables         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticDeformationGradientdPreviousStateVariables, floatVector, setdUpdatedPlasticDeformationGradientdPreviousStateVariables )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdMicroStress,               fourthOrderTensor, setdUpdatedPlasticMicroDeformationdMicroStress               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousMicroStress,       fourthOrderTensor, setdUpdatedPlasticMicroDeformationdPreviousMicroStress       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdF,                         fourthOrderTensor, setdUpdatedPlasticMicroDeformationdF                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousF,                 fourthOrderTensor, setdUpdatedPlasticMicroDeformationdPreviousF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdFn,                        floatVector, setdUpdatedPlasticMicroDeformationdFn                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousFn,                floatVector, setdUpdatedPlasticMicroDeformationdPreviousFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdChi,                       fourthOrderTensor, setdUpdatedPlasticMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousChi,               fourthOrderTensor, setdUpdatedPlasticMicroDeformationdPreviousChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdChin,                      floatVector, setdUpdatedPlasticMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousChin,              floatVector, setdUpdatedPlasticMicroDeformationdPreviousChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdStateVariables,            floatVector, setdUpdatedPlasticMicroDeformationdStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticMicroDeformationdPreviousStateVariables,    floatVector, setdUpdatedPlasticMicroDeformationdPreviousStateVariables    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdMacroStress,         fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdMacroStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress, fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousMacroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdMicroStress,         fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdMicroStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress, fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousMicroStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdHigherOrderStress,         sixthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdHigherOrderStress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress, sixthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdF,                         fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdF                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousF,                 fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdFn,                        floatVector, setdUpdatedPlasticGradientMicroDeformationdFn                        )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousFn,                floatVector, setdUpdatedPlasticGradientMicroDeformationdPreviousFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdChi,                       fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdChi                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousChi,               fifthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousChi               )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdChin,                      floatVector, setdUpdatedPlasticGradientMicroDeformationdChin                      )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousChin,              floatVector, setdUpdatedPlasticGradientMicroDeformationdPreviousChin              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdGradChi,                   sixthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdGradChi                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousGradChi,           sixthOrderTensor, setdUpdatedPlasticGradientMicroDeformationdPreviousGradChi           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdGradChin,                  floatVector, setdUpdatedPlasticGradientMicroDeformationdGradChin                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousGradChin,          floatVector, setdUpdatedPlasticGradientMicroDeformationdPreviousGradChin          )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdStateVariables,            floatVector, setdUpdatedPlasticGradientMicroDeformationdStateVariables            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables,    floatVector, setdUpdatedPlasticGradientMicroDeformationdPreviousStateVariables    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableResiduals,                                            floatVector, setStateVariableResiduals                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, stateVariableJacobians,                                            floatVector, setStateVariableJacobians                                            )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableResidualsdD,                                         floatVector, setdStateVariableResidualsdD                                         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dStateVariableResidualsdPreviousISVs,                              floatVector, setdStateVariableResidualsdPreviousISVs                              )

        };

    }

}

#endif
