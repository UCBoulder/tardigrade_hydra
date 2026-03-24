/**
 ******************************************************************************
 * \file tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.h
 ******************************************************************************
 * An implementation of micromorphic drucker-prager plasticity model which is
 * based on a radial return mapping approach to its solution.
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_MICROMORPHIC_RADIAL_RETURN_DRUCKER_PRAGER_PLASTICITY_H
#define TARDIGRADE_HYDRA_MICROMORPHIC_RADIAL_RETURN_DRUCKER_PRAGER_PLASTICITY_H

#define USE_EIGEN
#include <tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include <tardigrade_vector_tools.h>

namespace tardigradeHydra {

    namespace micromorphicRadialReturnDruckerPragerPlasticity {

        // forward class definitions
        namespace unit_test {
            class residualTester;
        }

        constexpr const char *str_end(const char *str) {
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
        constexpr const char *r_slant(const char *str) {
            /*! Recursively search string for rightmost UNIX path separator from the right
             * \param *str: pointer to string END of UNIX path like string
             * \return *str: pointer to start of base name
             */
            return *str == '/' ? (str + 1) : r_slant(str - 1);
        }
        constexpr const char *file_name(const char *str) {
            /*! Return the current file name with extension at compile time
             * \param *str: pointer to string START of UNIX path like string
             * \return str: file base name
             */
            return str_slant(str) ? r_slant(str_end(str)) : str;
        }
        // Return filename for constructing debugging messages
        // https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
        const std::string __BASENAME__ = file_name(__FILE__);  //!< The base filename which will be parsed
        const std::string __FILENAME__ =
            __BASENAME__.substr(0, __BASENAME__.find_last_of("."));  //!< The parsed filename for error handling

        typedef tardigradeErrorTools::Node           errorNode;    //!< Redefinition for the error node
        typedef errorNode                           *errorOut;     //!< Redefinition for a pointer to the error node
        typedef double                               floatType;    //!< Define the float values type.
        typedef std::vector<floatType>               floatVector;  //!< Define a vector of floats
        typedef std::vector<std::vector<floatType> > floatMatrix;  //!< Define a matrix of floats

        typedef floatType                               variableType;    //!< Define the variable values type.
        typedef std::vector<variableType>               variableVector;  //!< Define a vector of variables
        typedef std::vector<std::vector<variableType> > variableMatrix;  //!< Define a matrix of variables

        typedef double                     parameterType;    //!< Define the parameter values type.
        typedef std::vector<parameterType> parameterVector;  //!< Define a vector of parameters

        typedef double                                  constantType;    //!< Define the constant values type.
        typedef std::vector<constantType>               constantVector;  //!< Define a vector of constants
        typedef std::vector<std::vector<constantType> > constantMatrix;  //!< Define a matrix of constants

        typedef std::vector<floatType> seventhOrderTensor;  //!< A seventh order tensor
        typedef std::vector<floatType> eighthOrderTensor;   //!< A eighth order tensor

        /*!
         * The residual for a micromorphic Drucker Prager plasticity model
         */
        class residual : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual {
           public:
            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setResidual;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setJacobian;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setdRdD;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setdRdT;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStresses;

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStressesJacobians;

            residual(hydraBaseMicromorphic *_hydra, const unsigned int &_numEquations,
                     const unsigned int              &plasticConfigurationIndex,
                     const std::vector<unsigned int> &stateVariableIndices, const floatVector &parameters,
                     const floatType integrationParameter = 0.5, const double yieldTolerance = 1e-3,
                     const double plasticMultiplierTolerance = -1e-4)
                : tardigradeHydra::micromorphicDruckerPragerPlasticity::residual(_hydra, _numEquations,
                                                                                 plasticConfigurationIndex,
                                                                                 stateVariableIndices, parameters,
                                                                                 integrationParameter),
                  _yieldTolerance(yieldTolerance),
                  _plasticMultiplierTolerance(plasticMultiplierTolerance) {
                /*!
                 * The main initialization constructor for the Drucker Prager plasticity residual based on optimization
                 *
                 * \param *_hydra: A pointer to the containing hydra class
                 * \param &_numEquations: The number of equations the residual defines
                 * \param &plasticConfigurationIndex: The index of the configuration which represents the plastic
                 * deformation \param &stateVariableIndices: The indices of the plastic state variables \param
                 * &parameters: The parameter vector \param &integrationParameter: The integration parameter for the
                 * function. 0 is explicit, 1 is implicit. \param &yieldTolerance: The allowable value of the yield
                 * stress before it is assumed to be an inactive surface (defaults to 1e-3) \param
                 * &plasticMultiplierTolerance: The allowable value of a plastic multiplier before it is assumed to be
                 * an inactive surface (defaults to 1e-4)
                 */
            }

            virtual void successfulIterativeStep() override;

            virtual void preIterativeSolve() override;

            virtual void postIterativeSolve() override;

            virtual void setupRelaxedStep(const unsigned int &relaxedStep) override;

           protected:
            virtual void setActiveConstraints();

            virtual void setChangeInConstraints(const bool value) {
                /*!
                 * Set whether there has been a change in the constraints
                 *
                 * \param value: The incoming value
                 */

                _change_in_constraints = value;
            }

            virtual bool getChangeInConstraints() {
                /*!
                 * Get whether there has been a change in the constraints
                 */

                return _change_in_constraints;
            }

            std::vector<bool> *getMutableActiveConstraints() {
                /*!
                 * Return a mutable reference to the active constraints
                 */

                return &_activeConstraints.second;
            }

            virtual floatType getConstraintChangeToleranceScaleFactor() {
                /*!
                 * Get the scale factor for the global tolerance if the constraints have changed
                 */

                return 1.0;
            }

            virtual void updateActiveConstraints();

            virtual void projectSuggestedX(std::vector<floatType> &trialX, const std::vector<floatType> &Xp) override;

            virtual void setStateVariableResiduals() override;

            virtual void setStateVariableJacobians() override;

            virtual void setdStateVariableResidualsdD() override;

            virtual void correctResiduals();

            const double getYieldTolerance() {
                /*!
                 * Return the yield tolerance
                 */

                return _yieldTolerance;
            }

            const double getPlasticMultiplierTolerance() {
                /*!
                 * Return the plastic multiplier tolerance
                 */

                return _plasticMultiplierTolerance;
            }

           private:
            TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(private, activeConstraints, std::vector<bool>,
                                                      setActiveConstraints)

            //! The allowable tolerance on the yield function
            double _yieldTolerance;

            //! The allowable tolerance on the plastic multiplier
            double _plasticMultiplierTolerance;

            //! Whether the constraints changed in the previous iteration
            bool _change_in_constraints = false;
        };

    }  // namespace micromorphicRadialReturnDruckerPragerPlasticity

}  // namespace tardigradeHydra

#include "tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.tpp"

#endif
