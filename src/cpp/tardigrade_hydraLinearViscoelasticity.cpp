/**
 ******************************************************************************
 * \file tardigrade_hydraLinearViscoelasticity.cpp
 ******************************************************************************
 * An implementation of linear elasticity using the hydra framework. Used as an
 * example and as the basis for more complex models.
 ******************************************************************************
 */

#include <tardigrade_constitutive_tools.h>
#include <tardigrade_hydraLinearViscoelasticity.h>
#include <tardigrade_stress_tools.h>

namespace tardigradeHydra {

    namespace linearViscoelasticity {

        void residual::decomposeParameterVector(const floatVector &parameters) {
            /*!
             * Decompose the parameter vector
             *
             * \param &parameters: The parameter vector. Assumed to be a vector of at least length 10 which defines
             *   num_volumetric_terms, num_isochoric_terms, Kinf, Ginf, volumetric_temperature_terms,
             * isochoric_temperature_terms, Ks, Ktaus, Gs, Gtaus where num_volumetric_terms and num_isochoric_terms are
             * the number of volumetric and isochoric Maxwell elements respectively, Kinf and Ginf are the infinite bulk
             * and shear moduli, volumetric_temperature_terms and isochoric_temperature_terms are the volumetric
             * reference temperature, \f$ C_1 \f$ and \f$ C_2 \f$ terms for the WLF shift, Ks are the bulk moduli, Ktaus
             * are the volumetric time constants, Gs are the shear moduli, and Gtaus are the isochoric time constants,
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                parameters.size() >= 10,
                "Parameter vector is expected to have a length of at least 10 but has a length of " +
                    std::to_string(parameters.size()))

            setNumVolumetricViscousTerms((unsigned int)(parameters[0] + 0.5));

            setNumIsochoricViscousTerms((unsigned int)(parameters[1] + 0.5));

            setNumStateVariables(getNumVolumetricViscousTerms() + sot_dimension * getNumIsochoricViscousTerms());

            TARDIGRADE_ERROR_TOOLS_CHECK(
                getNumStateVariables() == (getViscoelasticISVUpperIndex() - getViscoelasticISVLowerIndex()),
                "The number of state variables required by the parameterization is not equal to the number of state "
                "variables indicated by the ISV bounds\n   required # ISVs: " +
                    std::to_string(getNumStateVariables()) +
                    "\n   ISV Lower Bound: " + std::to_string(getViscoelasticISVLowerIndex()) +
                    "\n   ISV Upper Bound: " + std::to_string(getViscoelasticISVUpperIndex()) + "\n")

            setKinf(parameters[2]);

            setGinf(parameters[3]);

            setVolumetricTemperatureParameters(floatVector(parameters.begin() + 4, parameters.begin() + 7));

            setIsochoricTemperatureParameters(floatVector(parameters.begin() + 7, parameters.begin() + 10));

            unsigned int parameterCount = 10 + 2 * getNumVolumetricViscousTerms() + 2 * getNumIsochoricViscousTerms();

            TARDIGRADE_ERROR_TOOLS_CHECK(parameters.size() == parameterCount,
                                         "The number of parameters provided is not consistent with the parameter "
                                         "counts\n  num parameters:      " +
                                             std::to_string(parameters.size()) + "\n  num viscous terms:   " +
                                             std::to_string(getNumVolumetricViscousTerms()) +
                                             "\n  num isochoric terms: " +
                                             std::to_string(getNumIsochoricViscousTerms()) +
                                             "\nThe number of parameters is 4 + 2 * ( numVolumetricViscousTerms + "
                                             "numIsochoricViscousTerms )\n  required parameter count: " +
                                             std::to_string(parameterCount) + "\n")

            unsigned int lb = 10;
            unsigned int ub = lb + getNumVolumetricViscousTerms();

            floatVector Ks(parameters.begin() + lb, parameters.begin() + ub);

            lb = ub;
            ub = lb + getNumVolumetricViscousTerms();

            floatVector Ktaus(parameters.begin() + lb, parameters.begin() + ub);

            lb = ub;
            ub = lb + getNumIsochoricViscousTerms();

            floatVector Gs(parameters.begin() + lb, parameters.begin() + ub);

            lb = ub;
            ub = lb + getNumIsochoricViscousTerms();

            floatVector Gtaus(parameters.begin() + lb, parameters.begin() + ub);

            setVolumetricModuli(Ks);

            setVolumetricTaus(Ktaus);

            setIsochoricModuli(Gs);

            setIsochoricTaus(Gtaus);
        }

        void residual::addParameterizationInfo(std::string &parameterization_info) {
            /*!
             * Add parameterization info to the incoming string
             *
             * \param &parameterization_info: The incoming string
             */

            std::stringstream ss;

            ss << "class: tardigradeHydra::linearViscoelasticity::residual\n\n";
            ss << "     name,                               description,       units, current value\n";
            ss << "    n_vol, The number of volumetric Maxwell elements,        none, "
               << getNumVolumetricViscousTerms() << "\n";
            ss << "    n_iso,  The number of isochoric Maxwell elements,        none, " << getNumIsochoricViscousTerms()
               << "\n";
            ss.precision(9);
            ss << std::scientific;
            ss << "     Kinf,           The infinite volumetric modulus,      stress, " << getKinf() << "\n";
            ss << "     Ginf,            The infinite isochoric modulus,      stress, " << getGinf() << "\n";
            ss << " Tref_vol,      The volumetric reference temperature, temperature, "
               << (*getVolumetricTemperatureParameters())[0] << "\n";
            ss << "   C1_vol,           The volumetric WLF C1 parameter,        none, "
               << (*getVolumetricTemperatureParameters())[1] << "\n";
            ss << "   C2_vol,           The volumetric WLF C2 parameter, temperature, "
               << (*getVolumetricTemperatureParameters())[2] << "\n";
            ss << " Tref_iso,       The isochoric reference temperature, temperature, "
               << (*getIsochoricTemperatureParameters())[0] << "\n";
            ss << "   C1_iso,            The isochoric WLF C1 parameter,        none, "
               << (*getIsochoricTemperatureParameters())[1] << "\n";
            ss << "   C2_iso,            The isochoric WLF C2 parameter, temperature, "
               << (*getIsochoricTemperatureParameters())[2] << "\n";
            ss << "\nVolumetric Maxwell elements tau (time), K (stress):\n";

            for (auto v = std::begin(*getVolumetricModuli()); v != std::end(*getVolumetricModuli()); ++v) {
                ss << (*getVolumetricTaus())[(unsigned int)(v - std::begin(*getVolumetricModuli()))] << ", " << *v
                   << "\n";
            }
            ss << "\nIsochoric Maxwell elements tau (time), G (stress):\n";

            for (auto v = std::begin(*getIsochoricModuli()); v != std::end(*getIsochoricModuli()); ++v) {
                ss << (*getIsochoricTaus())[(unsigned int)(v - std::begin(*getIsochoricModuli()))] << ", " << *v
                   << "\n";
            }

            ss.unsetf(std::ios_base::floatfield);

            parameterization_info.append(ss.str());
        }

        void residual::setNumVolumetricViscousTerms(const unsigned int &num) {
            /*!
             * Set the number of volumetric prony-series viscous terms
             *
             * \param &num: The number of volumetric Prony-series terms
             */

            _numVolumetricViscousTerms = num;
        }

        void residual::setNumIsochoricViscousTerms(const unsigned int &num) {
            /*!
             * Set the number of isochoric prony-series viscous terms
             *
             * \param &num: The number of isochorric Prony-series terms
             */

            _numIsochoricViscousTerms = num;
        }

        void residual::setKinf(const floatType &Kinf) {
            /*!
             * Set the infinite bulk modulus
             *
             * \param &Kinf: The infinite bulk modulus
             */

            _Kinf = Kinf;
        }

        void residual::setGinf(const floatType &Ginf) {
            /*!
             * Set the infinite shear modulus
             *
             * \param &Ginf: The infinite shear modulus
             */

            _Ginf = Ginf;
        }

        void residual::setVolumetricModuli(const floatVector &Ks) {
            /*!
             * Set the volumetric moduli
             *
             * \param &Ks: The bulk moduli
             */

            _Ks = Ks;
        }

        void residual::setIsochoricModuli(const floatVector &Gs) {
            /*!
             * Set the isochoric moduli
             *
             * \param &Gs: The isochoric moduli
             */

            _Gs = Gs;
        }

        void residual::setVolumetricTaus(const floatVector &taus) {
            /*!
             * Set the volumetric time constants
             *
             * \param &taus: The bulk time constants
             */

            _volumetricTaus = taus;
        }

        void residual::setIsochoricTaus(const floatVector &taus) {
            /*!
             * Set the isochoric time constants
             *
             * \param &taus: The isochoric time constants
             */

            _isochoricTaus = taus;
        }

        void residual::decomposeElasticDeformation() {
            /*!
             * Decompose the elastic deformation into volumetric and isochoric parts
             */

            auto Je = get_SetDataStorage_Je();

            auto Fehat = get_SetDataStorage_Fehat();

            TARDIGRADE_ERROR_TOOLS_CATCH(decomposeDeformation(*get_Fe(), *Je.value, *Fehat.value));
        }

        void residual::decomposePreviousElasticDeformation() {
            /*!
             * Decompose the previous elastic deformation into volumetric and isochoric parts
             */

            auto previousJe = get_SetDataStorage_previousJe();

            auto previousFehat = get_SetDataStorage_previousFehat();

            TARDIGRADE_ERROR_TOOLS_CATCH(decomposeDeformation(*get_previousFe(), *previousJe.value,
                                                              *previousFehat.value));
        }

        void residual::decomposeDeformation(const secondOrderTensor &F, floatType &J, secondOrderTensor &Fhat) {
            /*!
             * Decompose a deformation into volumetric and isochoric parts where
             *
             * \f$\hat{\bf{F}} = J^{-\frac{1}{3}} \bf{F} \f$
             *
             * \param &F: The incoming deformation gradient
             * \param &J: The jacobian of deformation
             * \param &Fhat: The isochoric part of the deformation gradient
             */

            using F_iterator = decltype(std::begin(F));
            using F_type     = std::iterator_traits<F_iterator>::value_type;

            TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(std::end(F) - std::begin(F)) == dimension * dimension,
                                         "The deformation has a size of " +
                                             std::to_string((unsigned int)(std::end(F) - std::begin(F))) +
                                             " but must have a size of " + std::to_string(dimension * dimension));

            J = tardigradeVectorTools::determinant<F_iterator, F_type, dimension, dimension>(std::begin(F), std::end(F),
                                                                                             dimension, dimension);

            TARDIGRADE_ERROR_TOOLS_CATCH(Fhat = F / std::pow(J, 1. / 3));
        }

        void residual::setdJedFe(const bool isPrevious) {
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             *
             * \param isPrevious: Flag for if the derivative is of the current (false) or previous (true) value
             */

            const secondOrderTensor *Fe;

            SetDataStorageBase<secondOrderTensor> dJedFe;

            if (isPrevious) {
                Fe = get_previousFe();

                dJedFe = get_SetDataStorage_previousdJedFe();

            } else {
                Fe = get_Fe();

                dJedFe = get_SetDataStorage_dJedFe();
            }

            using Fe_iterator     = decltype(std::begin(*Fe));
            using dJedFe_iterator = decltype(std::begin(*dJedFe.value));
            using Fe_type         = std::iterator_traits<Fe_iterator>::value_type;

            dJedFe.zero(dimension * dimension);

            tardigradeVectorTools::computeDDetADA<Fe_iterator, dJedFe_iterator, Fe_type, dimension, dimension>(
                std::begin(*Fe), std::end(*Fe), dimension, dimension, std::begin(*dJedFe.value),
                std::end(*dJedFe.value));
        }

        void residual::setdJedFe() {
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             */

            setdJedFe(false);
        }

        void residual::setdFehatdFe() {
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             */

            setdFehatdFe(false);
        }

        void residual::setPreviousdJedFe() {
            /*!
             * Set the previous derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             */

            setdJedFe(true);
        }

        void residual::setPreviousdFehatdFe() {
            /*!
             * Set the previous derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             */

            setdFehatdFe(true);
        }

        void residual::setdFehatdFe(const bool isPrevious) {
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             *
             * \param isPrevious: Flag for if the derivative is of the current (false) or previous (true) value
             */

            const floatType *Je;

            const secondOrderTensor *Fe;

            const secondOrderTensor *dJedFe;

            SetDataStorageBase<fourthOrderTensor> dFehatdFe;

            if (isPrevious) {
                Je = get_previousJe();

                Fe = get_previousFe();

                dJedFe = get_previousdJedFe();

                dFehatdFe = get_SetDataStorage_previousdFehatdFe();

            } else {
                Je = get_Je();

                Fe = get_Fe();

                dJedFe = get_dJedFe();

                dFehatdFe = get_SetDataStorage_dFehatdFe();
            }

            dFehatdFe.zero(sot_dimension * sot_dimension);
            for (unsigned int i = 0; i < sot_dimension; i++) {
                (*dFehatdFe.value)[sot_dimension * i + i] = std::pow((*Je), -1. / 3);
            }

            Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_Fe(Fe->data(), sot_dimension);

            Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_dJedFe(dJedFe->data(), sot_dimension);

            Eigen::Map<Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> > map_dFehatdFe(
                dFehatdFe.value->data(), sot_dimension, sot_dimension);

            map_dFehatdFe -= (map_Fe * map_dJedFe.transpose() * std::pow((*Je), -4. / 3) / 3).eval();
        }

        void residual::decomposeStateVariableVector(floatVector &volumetricStateVariables,
                                                    floatVector &isochoricStateVariables) {
            /*!
             * Decompose the state variable vector into parameters associated with the
             * volumetric and isochoric viscoelasticity
             */

            unsigned int lb = getViscoelasticISVLowerIndex();
            unsigned int ub = lb + getNumVolumetricViscousTerms();

            volumetricStateVariables = floatVector(hydra->get_additionalStateVariables()->begin() + lb,
                                                   hydra->get_additionalStateVariables()->begin() + ub);

            lb = ub;

            ub = lb + sot_dimension * getNumIsochoricViscousTerms();

            isochoricStateVariables = floatVector(hydra->get_additionalStateVariables()->begin() + lb,
                                                  hydra->get_additionalStateVariables()->begin() + ub);
        }

        void residual::updateAdditionalStateVariables(floatVector &additionalStateVariables) {
            /*!
             * Update the additional state variables with the values stored in the volumetric and
             * isochoric state variable arrays
             */

            std::copy(std::begin(*get_volumetricViscoelasticStateVariables()),
                      std::end(*get_volumetricViscoelasticStateVariables()),
                      std::begin(additionalStateVariables) + getViscoelasticISVLowerIndex());

            std::copy(std::begin(*get_isochoricViscoelasticStateVariables()),
                      std::end(*get_isochoricViscoelasticStateVariables()),
                      std::begin(additionalStateVariables) + getViscoelasticISVLowerIndex() +
                          getNumVolumetricViscousTerms());
        }

        void residual::setNumStateVariables(const unsigned int numStateVariables) {
            /*!
             * Set the number of state variables
             *
             * \param &numStateVariables: The number of state variables required
             */

            _numStateVariables = numStateVariables;
        }

        floatType residual::computeRateMultiplier(const floatVector &variables, const floatVector &parameters) {
            /*!
             * Compute the value of the rate multiplier this implementation uses the
             * WLF function of temperature to compute an increase in evolution rate for
             * high temperatures and a decrease for low temperatures.
             *
             * \param &variables: The incoming variables { temperature }
             * \param &parameters: The incoming parameters see tardigradeConstitutiveTools::WLF
             */

            if (variables.size() != 1) {
                TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("The incoming variables must have a size 1"));
            }

            floatType invRM;

            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::WLF(variables[0], parameters, invRM));

            return 1 / invRM;
        }

        floatVector residual::computedRateMultiplierdVariables(const floatVector &variables,
                                                               const floatVector &parameters) {
            /*!
             * Compute the value of the derivative of the rate multiplier with respect to
             * the variables. This implementation uses the WLF function of temperature to
             * compute an increase in evolution rate for high temperatures and a decrease
             * for low temperatures.
             *
             * \param &variables: The incoming variables { temperature }
             * \param &parameters: The incoming parameters see tardigradeConstitutiveTools::WLF
             */

            if (variables.size() != 1) {
                TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("The incoming variables must have a size 1"));
            }

            floatType   invRM;
            floatVector dinvRMdT(1, 0);

            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::WLF(variables[0], parameters, invRM,
                                                                          dinvRMdT[0]));

            return -1 / (invRM * invRM) * dinvRMdT;
        }

        void residual::setVolumetricRateMultiplier() {
            /*!
             * Set the value of the volumetric rate multiplier
             */

            auto rateMultiplier = get_SetDataStorage_volumetricRateMultiplier();

            TARDIGRADE_ERROR_TOOLS_CATCH(*rateMultiplier.value =
                                             computeRateMultiplier({hydra->getTemperature()},
                                                                   *getVolumetricTemperatureParameters()));
        }

        void residual::setPreviousVolumetricRateMultiplier() {
            /*!
             * Set the previous value of the volumetric rate multiplier
             */

            auto rateMultiplier = get_SetDataStorage_previousVolumetricRateMultiplier();

            TARDIGRADE_ERROR_TOOLS_CATCH(*rateMultiplier.value =
                                             computeRateMultiplier({hydra->getPreviousTemperature()},
                                                                   *getVolumetricTemperatureParameters()));
        }

        void residual::setIsochoricRateMultiplier() {
            /*!
             * Set the value of the isochoric rate multiplier
             */

            auto rateMultiplier = get_SetDataStorage_isochoricRateMultiplier();

            TARDIGRADE_ERROR_TOOLS_CATCH(*rateMultiplier.value =
                                             computeRateMultiplier({hydra->getTemperature()},
                                                                   *getIsochoricTemperatureParameters()));
        }

        void residual::setPreviousIsochoricRateMultiplier() {
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            auto rateMultiplier = get_SetDataStorage_previousIsochoricRateMultiplier();

            TARDIGRADE_ERROR_TOOLS_CATCH(*rateMultiplier.value =
                                             computeRateMultiplier({hydra->getPreviousTemperature()},
                                                                   *getIsochoricTemperatureParameters()));
        }

        void residual::setdVolumetricRateMultiplierdT() {
            /*!
             * Set the value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            auto dRateMultiplierdT = get_SetDataStorage_dVolumetricRateMultiplierdT();

            TARDIGRADE_ERROR_TOOLS_CATCH(*dRateMultiplierdT.value = computedRateMultiplierdVariables(
                                             {hydra->getTemperature()}, *getVolumetricTemperatureParameters())[0]);
        }

        void residual::setdPreviousVolumetricRateMultiplierdPreviousT() {
            /*!
             * Set the previous value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            auto dRateMultiplierdT = get_SetDataStorage_dPreviousVolumetricRateMultiplierdPreviousT();

            TARDIGRADE_ERROR_TOOLS_CATCH(*dRateMultiplierdT.value = computedRateMultiplierdVariables(
                                             {hydra->getPreviousTemperature()},
                                             *getVolumetricTemperatureParameters())[0]);
        }

        void residual::setdIsochoricRateMultiplierdT() {
            /*!
             * Set the value of the derivative of the isochoric rate multiplier
             * with respect to the temperature
             */

            auto dRateMultiplierdT = get_SetDataStorage_dIsochoricRateMultiplierdT();

            TARDIGRADE_ERROR_TOOLS_CATCH(*dRateMultiplierdT.value =
                                             computedRateMultiplierdVariables({hydra->getTemperature()},
                                                                              *getIsochoricTemperatureParameters())[0]);
        }

        void residual::setdPreviousIsochoricRateMultiplierdPreviousT() {
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            auto dRateMultiplierdT = get_SetDataStorage_dPreviousIsochoricRateMultiplierdPreviousT();

            TARDIGRADE_ERROR_TOOLS_CATCH(*dRateMultiplierdT.value =
                                             computedRateMultiplierdVariables({hydra->getPreviousTemperature()},
                                                                              *getIsochoricTemperatureParameters())[0]);
        }

        void residual::setVolumetricTemperatureParameters(const floatVector &parameters) {
            /*!
             * Set the volumetric temperature parameters
             *
             * \param &parameters: The parameters for the volumetric temperature dependence
             */

            _volumetricTemperatureParameters = parameters;
        }

        void residual::setIsochoricTemperatureParameters(const floatVector &parameters) {
            /*!
             * Set the isochoric temperature parameters
             *
             * \param &parameters: The parameters for the isochoric temperature dependence
             */

            _isochoricTemperatureParameters = parameters;
        }

        floatVector residual::getVolumetricViscoelasticParameters() {
            /*!
             * Get the volumetric viscoelastic parameters prepared for tardigradeStressTools::linearViscoelasticity
             */

            floatVector parameters(1 + 2 * getNumVolumetricViscousTerms(), 0);

            parameters[0] = getKinf();

            for (unsigned int i = 1; i <= getNumVolumetricViscousTerms(); i++) {
                parameters[i]                                  = (*getVolumetricTaus())[i - 1];
                parameters[i + getNumVolumetricViscousTerms()] = (*getVolumetricModuli())[i - 1];
            }

            return parameters;
        }

        floatVector residual::getIsochoricViscoelasticParameters() {
            /*!
             * Get the isochoric viscoelastic parameters prepared for tardigradeStressTools::linearViscoelasticity
             */

            floatVector parameters(1 + 2 * getNumIsochoricViscousTerms(), 0);

            parameters[0] = 2 * getGinf();

            for (unsigned int i = 1; i <= getNumIsochoricViscousTerms(); i++) {
                parameters[i]                                 = (*getIsochoricTaus())[i - 1];
                parameters[i + getNumIsochoricViscousTerms()] = 2 * (*getIsochoricModuli())[i - 1];
            }

            return parameters;
        }

        void residual::setPK2MeanStress() {
            /*!
             * Set the mean stress for the Second Piola-Kirchhoff stress
             */

            setPK2MeanStress(false);
        }

        void residual::setPreviousPK2MeanStress() {
            /*!
             * Set the previous mean stress for the Second Piola-Kirchhoff stress
             */

            setPK2MeanStress(true);
        }

        void residual::setPK2MeanStress(const bool isPrevious) {
            /*!
             * Set the mean stress for the Second Piola-Kirchhoff stress
             *
             * \param &isPrevious: Flag for if the previous (true) or current (false) stress should be calculated
             */

            const floatType *Je;

            const secondOrderTensor *dJedFe;

            const floatType *previousJe = get_previousJe();

            floatType time;

            floatType previousTime = hydra->getTime() - hydra->getDeltaTime();

            const floatType *volumetricRateMultiplier;

            const floatType *dVolumetricRateMultiplierdT;

            const floatType *previousVolumetricRateMultiplier = get_previousVolumetricRateMultiplier();

            SetDataStorageBase<floatType> PK2MeanStress;

            SetDataStorageBase<secondOrderTensor> dPK2MeanStressdFe;

            SetDataStorageBase<floatType> dPK2MeanStressdT;

            if (isPrevious) {
                time = previousTime;

                Je = get_previousJe();

                dJedFe = get_previousdJedFe();

                volumetricRateMultiplier = get_previousVolumetricRateMultiplier();

                dVolumetricRateMultiplierdT = get_dPreviousVolumetricRateMultiplierdPreviousT();

                PK2MeanStress = get_SetDataStorage_previousPK2MeanStress();

                dPK2MeanStressdFe = get_SetDataStorage_previousdPK2MeanStressdFe();

                dPK2MeanStressdT = get_SetDataStorage_previousdPK2MeanStressdT();

            } else {
                time = hydra->getTime();

                Je = get_Je();

                dJedFe = get_dJedFe();

                volumetricRateMultiplier = get_volumetricRateMultiplier();

                dVolumetricRateMultiplierdT = get_dVolumetricRateMultiplierdT();

                PK2MeanStress = get_SetDataStorage_PK2MeanStress();

                dPK2MeanStressdFe = get_SetDataStorage_dPK2MeanStressdFe();

                dPK2MeanStressdT = get_SetDataStorage_dPK2MeanStressdT();
            }

            // Compute the strain measures

            floatVector volumetricStrain = {(*Je - 1)};

            floatVector previousVolumetricStrain = {(*previousJe - 1)};

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH(decomposeStateVariableVector(previousVolumetricStateVariables,
                                                                      previousIsochoricStateVariables));

            floatVector _PK2MeanStress;

            floatVector deltaPK2MeanStress;

            floatMatrix _dPK2MeanStressdJe;

            floatVector dPK2MeanStressdJe;

            floatVector dPK2MeanStressdRateModifier;

            floatMatrix _dPK2MeanStressdPreviousJe;

            floatVector dPK2MeanStressdPreviousJe;

            floatVector dPK2MeanStressdPreviousRateModifier;

            floatMatrix _dPK2MeanStressdPreviousVolumetricISVs;

            floatVector dPK2MeanStressdPreviousVolumetricISVs;

            floatMatrix _dISVsdJe;

            floatVector dISVsdJe;

            floatVector dISVsdRateModifier;

            floatMatrix _dISVsdPreviousJe;

            floatVector dISVsdPreviousJe;

            floatVector dISVsdPreviousRateModifier;

            floatMatrix _dISVsdPreviousVolumetricISVs;

            floatVector dISVsdPreviousVolumetricISVs;

            floatVector currentVolumetricStateVariables;

            // Compute the viscous mean stress

            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeStressTools::linearViscoelasticity(
                time, volumetricStrain, previousTime, previousVolumetricStrain, *volumetricRateMultiplier,
                *previousVolumetricRateMultiplier, previousVolumetricStateVariables,
                getVolumetricViscoelasticParameters(), getIntegrationAlpha(), deltaPK2MeanStress, _PK2MeanStress,
                currentVolumetricStateVariables, _dPK2MeanStressdJe, dPK2MeanStressdRateModifier,
                _dPK2MeanStressdPreviousJe, dPK2MeanStressdPreviousRateModifier, _dPK2MeanStressdPreviousVolumetricISVs,
                _dISVsdJe, dISVsdRateModifier, _dISVsdPreviousJe, dISVsdPreviousRateModifier,
                _dISVsdPreviousVolumetricISVs));

            *PK2MeanStress.value = _PK2MeanStress[0];

            *dPK2MeanStressdFe.value = _dPK2MeanStressdJe[0][0] * (*dJedFe);

            *dPK2MeanStressdT.value = dPK2MeanStressdRateModifier[0] * (*dVolumetricRateMultiplierdT);

            dPK2MeanStressdPreviousJe = tardigradeVectorTools::appendVectors(_dPK2MeanStressdPreviousJe);

            dPK2MeanStressdPreviousVolumetricISVs =
                tardigradeVectorTools::appendVectors(_dPK2MeanStressdPreviousVolumetricISVs);

            dISVsdJe = tardigradeVectorTools::appendVectors(_dISVsdJe);

            dISVsdPreviousJe = tardigradeVectorTools::appendVectors(_dISVsdPreviousJe);

            dISVsdPreviousVolumetricISVs = tardigradeVectorTools::appendVectors(_dISVsdPreviousVolumetricISVs);

            if (!isPrevious) {
                set_volumetricViscoelasticStateVariables(currentVolumetricStateVariables);

                set_dPK2MeanStressdPreviousFe(dPK2MeanStressdPreviousJe[0] * (*get_previousdJedFe()));

                set_dPK2MeanStressdPreviousT(dPK2MeanStressdPreviousRateModifier[0] *
                                             (*get_dPreviousVolumetricRateMultiplierdPreviousT()));

                const unsigned int vol_isvs_size = previousVolumetricStateVariables.size();

                const unsigned int iso_isvs_size = previousIsochoricStateVariables.size();

                auto dPK2MeanStressdPreviousISVs = get_SetDataStorage_dPK2MeanStressdPreviousISVs();
                dPK2MeanStressdPreviousISVs.zero(vol_isvs_size + iso_isvs_size);

                auto dVolumetricISVsdPreviousISVs = get_SetDataStorage_dVolumetricISVsdPreviousISVs();
                dVolumetricISVsdPreviousISVs.zero(vol_isvs_size * (vol_isvs_size + iso_isvs_size));

                for (unsigned int i = 0; i < vol_isvs_size; i++) {
                    (*dPK2MeanStressdPreviousISVs.value)[i] = dPK2MeanStressdPreviousVolumetricISVs[i];

                    for (unsigned int j = 0; j < vol_isvs_size; j++) {
                        (*dVolumetricISVsdPreviousISVs.value)[(vol_isvs_size + iso_isvs_size) * i + j] =
                            dISVsdPreviousVolumetricISVs[(vol_isvs_size)*i + j];
                    }
                }

                Eigen::Map<const Eigen::Vector<floatType, -1> > map_dISVsdJe(dISVsdJe.data(), dISVsdJe.size());
                Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_dJedFe(dJedFe->data(), sot_dimension);

                auto dVolumetricISVsdFe = get_SetDataStorage_dVolumetricISVsdFe();
                dVolumetricISVsdFe.zero(dISVsdJe.size() * sot_dimension);
                Eigen::Map<Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> > map_dVolumetricISVsdFe(
                    dVolumetricISVsdFe.value->data(), dISVsdJe.size(), sot_dimension);

                map_dVolumetricISVsdFe = (map_dISVsdJe * map_dJedFe.transpose()).eval();

                auto dVolumetricISVsdT   = get_SetDataStorage_dVolumetricISVsdT();
                *dVolumetricISVsdT.value = dISVsdRateModifier * (*dVolumetricRateMultiplierdT);

                Eigen::Map<const Eigen::Vector<floatType, -1> > map_dISVsdPreviousJe(dISVsdPreviousJe.data(),
                                                                                     dISVsdPreviousJe.size());
                Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_previousdJedFe(
                    get_previousdJedFe()->data(), sot_dimension);

                auto dVolumetricISVsdPreviousFe = get_SetDataStorage_dVolumetricISVsdPreviousFe();
                dVolumetricISVsdPreviousFe.zero(dISVsdPreviousJe.size() * sot_dimension);
                Eigen::Map<Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> >
                    map_dVolumetricISVsdPreviousFe(dVolumetricISVsdPreviousFe.value->data(), dISVsdPreviousJe.size(),
                                                   sot_dimension);

                map_dVolumetricISVsdPreviousFe = (map_dISVsdPreviousJe * map_previousdJedFe.transpose()).eval();

                auto dVolumetricISVsdPreviousT = get_SetDataStorage_dVolumetricISVsdPreviousT();
                *dVolumetricISVsdPreviousT.value =
                    dISVsdPreviousRateModifier * (*get_dPreviousVolumetricRateMultiplierdPreviousT());
            }
        }

        void residual::setdPK2MeanStressdT() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the temperature
             */

            get_PK2MeanStress();
        }

        void residual::setdPK2MeanStressdFe() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation gradient
             */

            get_PK2MeanStress();
        }

        void residual::setdPK2MeanStressdPreviousT() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous temperature
             */

            get_PK2MeanStress();
        }

        void residual::setdPK2MeanStressdPreviousFe() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous elastic deformation gradient
             */

            get_PK2MeanStress();
        }

        void residual::setdPK2MeanStressdPreviousISVs() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous state variables
             */

            get_PK2MeanStress();
        }

        void residual::setdVolumetricISVsdT() {
            /*!
             * Set the derivative of the volumetric isvs w.r.t. the temperature
             */

            get_PK2MeanStress();
        }

        void residual::setdVolumetricISVsdFe() {
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the elastic deformation gradient
             */

            get_PK2MeanStress();
        }

        void residual::setdVolumetricISVsdPreviousT() {
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous temperature
             */

            get_PK2MeanStress();
        }

        void residual::setdVolumetricISVsdPreviousFe() {
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous elastic deformation gradient
             */

            get_PK2MeanStress();
        }

        void residual::setdVolumetricISVsdPreviousISVs() {
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous state variables
             */

            get_PK2MeanStress();
        }

        void residual::setPreviousdPK2MeanStressdT() {
            /*!
             * Set the previous derivative of the PK2 mean stress w.r.t. the temperature
             */

            get_previousPK2MeanStress();
        }

        void residual::setPreviousdPK2MeanStressdFe() {
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation gradient
             */

            get_previousPK2MeanStress();
        }

        void residual::setPK2IsochoricStress() {
            /*!
             * Set the isochoric second Piola-Kirchhoff stress
             */

            setPK2IsochoricStress(false);
        }

        void residual::setPreviousPK2IsochoricStress() {
            /*!
             * Set the previous isochoric second Piola-Kirchhoff stress
             */

            setPK2IsochoricStress(true);
        }

        void residual::setPK2IsochoricStress(const bool isPrevious) {
            /*!
             * Set the isochoric second Piola-Kirchhoff stress
             *
             * \param &isPrevious: Flag for if the previous (true) or current (false) stress should be calculated
             */

            const secondOrderTensor *Fehat;

            const secondOrderTensor *previousFehat = get_previousFehat();

            const fourthOrderTensor *dFehatdFe;

            floatType time;

            floatType previousTime = hydra->getTime() - hydra->getDeltaTime();

            const floatType *isochoricRateMultiplier;

            const floatType *dIsochoricRateMultiplierdT;

            const floatType *previousIsochoricRateMultiplier = get_previousIsochoricRateMultiplier();

            SetDataStorageBase<secondOrderTensor> PK2IsochoricStress;

            SetDataStorageBase<fourthOrderTensor> dPK2IsochoricStressdFe;

            SetDataStorageBase<secondOrderTensor> dPK2IsochoricStressdT;

            if (isPrevious) {
                time = previousTime;

                Fehat = get_previousFehat();

                dFehatdFe = get_previousdFehatdFe();

                isochoricRateMultiplier = get_previousIsochoricRateMultiplier();

                dIsochoricRateMultiplierdT = get_dPreviousIsochoricRateMultiplierdPreviousT();

                PK2IsochoricStress = get_SetDataStorage_previousPK2IsochoricStress();

                dPK2IsochoricStressdFe = get_SetDataStorage_previousdPK2IsochoricStressdFe();

                dPK2IsochoricStressdT = get_SetDataStorage_previousdPK2IsochoricStressdT();

            } else {
                time = hydra->getTime();

                Fehat = get_Fehat();

                dFehatdFe = get_dFehatdFe();

                isochoricRateMultiplier = get_isochoricRateMultiplier();

                dIsochoricRateMultiplierdT = get_dIsochoricRateMultiplierdT();

                PK2IsochoricStress = get_SetDataStorage_PK2IsochoricStress();

                dPK2IsochoricStressdFe = get_SetDataStorage_dPK2IsochoricStressdFe();

                dPK2IsochoricStressdT = get_SetDataStorage_dPK2IsochoricStressdT();
            }

            // Compute the strain measures

            secondOrderTensor isochoricStrain, previousIsochoricStrain;
            fourthOrderTensor dEehatdFehat, previousdEehatdFehat;
            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::computeGreenLagrangeStrain(*Fehat,
                                                                                                 isochoricStrain,
                                                                                                 dEehatdFehat));

            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::computeGreenLagrangeStrain(
                *previousFehat, previousIsochoricStrain, previousdEehatdFehat));

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> > map_dEehatdFehat(
                dEehatdFehat.data(), sot_dimension, sot_dimension);
            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> > map_dFehatdFe(
                dFehatdFe->data(), sot_dimension, sot_dimension);

            fourthOrderTensor dEehatdFe(sot_dimension * sot_dimension, 0);
            Eigen::Map<Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> > map_dEehatdFe(
                dEehatdFe.data(), sot_dimension, sot_dimension);

            map_dEehatdFe = (map_dEehatdFehat * map_dFehatdFe).eval();

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH(decomposeStateVariableVector(previousVolumetricStateVariables,
                                                                      previousIsochoricStateVariables));

            secondOrderTensor deltaPK2IsochoricStress;

            floatVector currentIsochoricStateVariables;

            floatMatrix _dPK2IsochoricStressdEe;

            fourthOrderTensor dPK2IsochoricStressdEe;

            secondOrderTensor dPK2IsochoricStressdRateMultiplier;

            floatMatrix _dPK2IsochoricStressdPreviousEe;

            fourthOrderTensor dPK2IsochoricStressdPreviousEe;

            secondOrderTensor dPK2IsochoricStressdPreviousRateMultiplier;

            floatMatrix _dPK2IsochoricStressdPreviousIsochoricISVs;

            floatVector dPK2IsochoricStressdPreviousIsochoricISVs;

            floatMatrix _dISVsdEe;

            floatVector dISVsdEe;

            floatVector dISVsdRateMultiplier;

            floatMatrix _dISVsdPreviousEe;

            floatVector dISVsdPreviousEe;

            floatVector dISVsdPreviousRateMultiplier;

            floatMatrix _dISVsdPreviousIsochoricISVs;

            floatVector dISVsdPreviousIsochoricISVs;

            // Compute the viscous isochoric stress
            TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeStressTools::linearViscoelasticity(
                time, isochoricStrain, previousTime, previousIsochoricStrain, *isochoricRateMultiplier,
                *previousIsochoricRateMultiplier, previousIsochoricStateVariables, getIsochoricViscoelasticParameters(),
                getIntegrationAlpha(), deltaPK2IsochoricStress, *PK2IsochoricStress.value,
                currentIsochoricStateVariables, _dPK2IsochoricStressdEe, dPK2IsochoricStressdRateMultiplier,
                _dPK2IsochoricStressdPreviousEe, dPK2IsochoricStressdPreviousRateMultiplier,
                _dPK2IsochoricStressdPreviousIsochoricISVs, _dISVsdEe, dISVsdRateMultiplier, _dISVsdPreviousEe,
                dISVsdPreviousRateMultiplier, _dISVsdPreviousIsochoricISVs));

            dPK2IsochoricStressdEe = tardigradeVectorTools::appendVectors(_dPK2IsochoricStressdEe);

            dPK2IsochoricStressdPreviousEe = tardigradeVectorTools::appendVectors(_dPK2IsochoricStressdPreviousEe);

            dPK2IsochoricStressdPreviousIsochoricISVs =
                tardigradeVectorTools::appendVectors(_dPK2IsochoricStressdPreviousIsochoricISVs);

            dISVsdEe = tardigradeVectorTools::appendVectors(_dISVsdEe);

            dISVsdPreviousEe = tardigradeVectorTools::appendVectors(_dISVsdPreviousEe);

            dISVsdPreviousIsochoricISVs = tardigradeVectorTools::appendVectors(_dISVsdPreviousIsochoricISVs);

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dPK2IsochoricStressdEe(dPK2IsochoricStressdEe.data(), sot_dimension, sot_dimension);

            dPK2IsochoricStressdFe.zero(sot_dimension * sot_dimension);
            Eigen::Map<Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dPK2IsochoricStressdFe(dPK2IsochoricStressdFe.value->data(), sot_dimension, sot_dimension);

            map_dPK2IsochoricStressdFe = (map_dPK2IsochoricStressdEe * map_dEehatdFe).eval();

            *dPK2IsochoricStressdT.value = dPK2IsochoricStressdRateMultiplier * (*dIsochoricRateMultiplierdT);

            if (!isPrevious) {
                set_isochoricViscoelasticStateVariables(currentIsochoricStateVariables);

                Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                    map_previousdEehatdFehat(previousdEehatdFehat.data(), sot_dimension, sot_dimension);

                Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                    map_previousdFehatdFe(get_previousdFehatdFe()->data(), sot_dimension, sot_dimension);

                fourthOrderTensor previousdEehatdFe(sot_dimension * sot_dimension, 0);

                Eigen::Map<Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                    map_previousdEehatdFe(previousdEehatdFe.data(), sot_dimension, sot_dimension);

                map_previousdEehatdFe = (map_previousdEehatdFehat * map_previousdFehatdFe).eval();

                auto dPK2IsochoricStressdPreviousFe = get_SetDataStorage_dPK2IsochoricStressdPreviousFe();

                dPK2IsochoricStressdPreviousFe.zero(sot_dimension * sot_dimension);

                Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                    map_dPK2IsochoricStressdPreviousEe(dPK2IsochoricStressdPreviousEe.data(), sot_dimension,
                                                       sot_dimension);

                Eigen::Map<Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                    map_dPK2IsochoricStressdPreviousFe(dPK2IsochoricStressdPreviousFe.value->data(), sot_dimension,
                                                       sot_dimension);

                map_dPK2IsochoricStressdPreviousFe =
                    (map_dPK2IsochoricStressdPreviousEe * map_previousdEehatdFe).eval();

                auto dPK2IsochoricStressdPreviousT = get_SetDataStorage_dPK2IsochoricStressdPreviousT();

                *dPK2IsochoricStressdPreviousT.value =
                    dPK2IsochoricStressdPreviousRateMultiplier * (*get_dPreviousIsochoricRateMultiplierdPreviousT());

                const unsigned int vol_isvs_size = previousVolumetricStateVariables.size();

                const unsigned int iso_isvs_size = previousIsochoricStateVariables.size();

                auto dPK2IsochoricStressdPreviousISVs = get_SetDataStorage_dPK2IsochoricStressdPreviousISVs();
                dPK2IsochoricStressdPreviousISVs.zero(sot_dimension * (vol_isvs_size + iso_isvs_size));

                auto dIsochoricISVsdPreviousISVs = get_SetDataStorage_dIsochoricISVsdPreviousISVs();
                dIsochoricISVsdPreviousISVs.zero(iso_isvs_size * (vol_isvs_size + iso_isvs_size));

                for (unsigned int i = 0; i < sot_dimension; i++) {
                    for (unsigned int j = 0; j < iso_isvs_size; j++) {
                        (*dPK2IsochoricStressdPreviousISVs
                              .value)[(vol_isvs_size + iso_isvs_size) * i + j + vol_isvs_size] =
                            dPK2IsochoricStressdPreviousIsochoricISVs[iso_isvs_size * i + j];
                    }
                }

                for (unsigned int i = 0; i < iso_isvs_size; i++) {
                    for (unsigned int j = 0; j < iso_isvs_size; j++) {
                        (*dIsochoricISVsdPreviousISVs.value)[(vol_isvs_size + iso_isvs_size) * i + j + vol_isvs_size] =
                            dISVsdPreviousIsochoricISVs[iso_isvs_size * i + j];
                    }
                }

                Eigen::Map<const Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> > map_dISVsdEe(
                    dISVsdEe.data(), iso_isvs_size, sot_dimension);

                auto dIsochoricISVsdFe = get_SetDataStorage_dIsochoricISVsdFe();
                dIsochoricISVsdFe.zero(iso_isvs_size * sot_dimension);
                Eigen::Map<Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> > map_dIsochoricISVsdFe(
                    dIsochoricISVsdFe.value->data(), iso_isvs_size, sot_dimension);

                map_dIsochoricISVsdFe = (map_dISVsdEe * map_dEehatdFe).eval();

                auto dIsochoricISVsdT   = get_SetDataStorage_dIsochoricISVsdT();
                *dIsochoricISVsdT.value = dISVsdRateMultiplier * (*dIsochoricRateMultiplierdT);

                Eigen::Map<const Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> > map_dISVsdPreviousEe(
                    dISVsdPreviousEe.data(), iso_isvs_size, sot_dimension);

                auto dIsochoricISVsdPreviousFe = get_SetDataStorage_dIsochoricISVsdPreviousFe();
                dIsochoricISVsdPreviousFe.zero(iso_isvs_size * sot_dimension);
                Eigen::Map<Eigen::Matrix<floatType, -1, sot_dimension, Eigen::RowMajor> > map_dIsochoricISVsdPreviousFe(
                    dIsochoricISVsdPreviousFe.value->data(), iso_isvs_size, sot_dimension);

                map_dIsochoricISVsdPreviousFe = (map_dISVsdPreviousEe * map_previousdEehatdFe).eval();

                auto dIsochoricISVsdPreviousT = get_SetDataStorage_dIsochoricISVsdPreviousT();
                *dIsochoricISVsdPreviousT.value =
                    dISVsdPreviousRateMultiplier * (*get_dPreviousIsochoricRateMultiplierdPreviousT());
            }
        }

        void residual::setdPK2IsochoricStressdFe() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            get_PK2IsochoricStress();
        }

        void residual::setdPK2IsochoricStressdT() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            get_PK2IsochoricStress();
        }

        void residual::setdPK2IsochoricStressdPreviousFe() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous elastic
             * deformation gradient
             */

            get_PK2IsochoricStress();
        }

        void residual::setdPK2IsochoricStressdPreviousT() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous temperature
             */

            get_PK2IsochoricStress();
        }

        void residual::setdPK2IsochoricStressdPreviousISVs() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous state variables
             */

            get_PK2IsochoricStress();
        }

        void residual::setdIsochoricISVsdFe() {
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the elastic
             * deformation gradient
             */

            get_PK2IsochoricStress();
        }

        void residual::setdIsochoricISVsdT() {
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the temperature
             */

            get_PK2IsochoricStress();
        }

        void residual::setdIsochoricISVsdPreviousFe() {
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous elastic
             * deformation gradient
             */

            get_PK2IsochoricStress();
        }

        void residual::setdIsochoricISVsdPreviousT() {
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous temperature
             */

            get_PK2IsochoricStress();
        }

        void residual::setdIsochoricISVsdPreviousISVs() {
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous state variables
             */

            get_PK2IsochoricStress();
        }

        void residual::setPreviousdPK2IsochoricStressdFe() {
            /*!
             * Set the previous derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            get_previousPK2IsochoricStress();
        }

        void residual::setPreviousdPK2IsochoricStressdT() {
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            get_previousPK2IsochoricStress();
        }

        void residual::setUpdatedVolumetricViscoelasticStateVariables() {
            /*!
             * Set the updated values of the volumetric viscoelastic state variables
             */

            get_PK2MeanStress();
        }

        void residual::setUpdatedIsochoricViscoelasticStateVariables() {
            /*!
             * Set the updated values of the isochoric viscoelastic state variables
             */

            get_PK2IsochoricStress();
        }

        void residual::setCurrentAdditionalStateVariables() {
            /*!
             * Set the updated current additional state variables
             */

            auto currentAdditionalStateVariables = get_SetDataStorage_currentAdditionalStateVariables();

            TARDIGRADE_ERROR_TOOLS_CATCH(*currentAdditionalStateVariables.value = tardigradeVectorTools::appendVectors(
                                             {*get_volumetricViscoelasticStateVariables(),
                                              *get_isochoricViscoelasticStateVariables()}));
        }

        void residual::setPK2Stress(const bool isPrevious) {
            /*!
             * Set the PK2 stress
             *
             * \param isPrevious: Flag for if to compute the current (false) or previous (true) PK2 stress
             */

            const secondOrderTensor *isochoric;

            const floatType *mean;

            SetDataStorageBase<secondOrderTensor> PK2Stress;

            if (isPrevious) {
                isochoric = get_previousPK2IsochoricStress();

                mean = get_previousPK2MeanStress();

                PK2Stress = get_SetDataStorage_previousPK2Stress();

            } else {
                isochoric = get_PK2IsochoricStress();

                mean = get_PK2MeanStress();

                PK2Stress = get_SetDataStorage_PK2Stress();
            }

            *PK2Stress.value = *isochoric;

            for (unsigned int i = 0; i < dimension; i++) {
                (*PK2Stress.value)[dimension * i + i] += *mean;
            }
        }

        void residual::setdPK2StressdFe() {
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            setdPK2StressdFe(false);
        }

        void residual::setPreviousdPK2StressdFe() {
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            setdPK2StressdFe(true);
        }

        void residual::setdPK2StressdFe(const bool isPrevious) {
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             *
             * \param isPrevious: Flag for whether to compute the derivative of the current (false) or previous (true)
             * stress
             */

            const fourthOrderTensor *dIsodFe;

            const secondOrderTensor *dMeandFe;

            SetDataStorageBase<fourthOrderTensor> dPK2StressdFe;

            if (isPrevious) {
                dIsodFe = get_previousdPK2IsochoricStressdFe();

                dMeandFe = get_previousdPK2MeanStressdFe();

                dPK2StressdFe = get_SetDataStorage_previousdPK2StressdFe();

            } else {
                dIsodFe = get_dPK2IsochoricStressdFe();

                dMeandFe = get_dPK2MeanStressdFe();

                dPK2StressdFe = get_SetDataStorage_dPK2StressdFe();
            }

            *dPK2StressdFe.value = *dIsodFe;

            for (unsigned int i = 0; i < dimension; i++) {
                for (unsigned int jk = 0; jk < sot_dimension; jk++) {
                    (*dPK2StressdFe.value)[dimension * dimension * dimension * i + dimension * dimension * i + jk] +=
                        (*dMeandFe)[jk];
                }
            }

            if (!isPrevious) {
                auto dPK2StressdPreviousFe = get_SetDataStorage_dPK2StressdPreviousFe();

                *dPK2StressdPreviousFe.value = *get_dPK2IsochoricStressdPreviousFe();

                for (unsigned int i = 0; i < dimension; i++) {
                    for (unsigned int jk = 0; jk < sot_dimension; jk++) {
                        (*dPK2StressdPreviousFe
                              .value)[dimension * dimension * dimension * i + dimension * dimension * i + jk] +=
                            (*get_dPK2MeanStressdPreviousFe())[jk];
                    }
                }
            }
        }

        void residual::setdPK2StressdPreviousFe() {
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous elastic deformation gradient
             */

            setdPK2StressdFe(false);
        }

        void residual::setdPK2StressdT() {
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            auto dPK2StressdT = get_SetDataStorage_dPK2StressdT();

            *dPK2StressdT.value = *get_dPK2IsochoricStressdT();

            for (unsigned int i = 0; i < dimension; i++) {
                (*dPK2StressdT.value)[dimension * i + i] += (*get_dPK2MeanStressdT());
            }
        }

        void residual::setdPK2StressdPreviousT() {
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous temperature
             */

            auto dPK2StressdPreviousT = get_SetDataStorage_dPK2StressdPreviousT();

            *dPK2StressdPreviousT.value = *get_dPK2IsochoricStressdPreviousT();

            for (unsigned int i = 0; i < dimension; i++) {
                (*dPK2StressdPreviousT.value)[dimension * i + i] += (*get_dPK2MeanStressdPreviousT());
            }
        }

        void residual::setPreviousdPK2StressdT() {
            /*!
             * Set the prevoius derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            auto previousdPK2StressdT = get_SetDataStorage_previousdPK2StressdT();

            *previousdPK2StressdT.value = *get_previousdPK2IsochoricStressdT();

            for (unsigned int i = 0; i < dimension; i++) {
                (*previousdPK2StressdT.value)[dimension * i + i] += (*get_previousdPK2MeanStressdT());
            }
        }

        void residual::setdPK2StressdPreviousISVs() {
            /*!
             * Set the prevoius derivative of the second Piola-Kirchhoff stress w.r.t. the previous ISVs
             */

            const unsigned int num_isvs = get_dPK2MeanStressdPreviousISVs()->size();

            auto dPK2StressdPreviousISVs = get_SetDataStorage_dPK2StressdPreviousISVs();

            *dPK2StressdPreviousISVs.value = *get_dPK2IsochoricStressdPreviousISVs();

            for (unsigned int i = 0; i < dimension; i++) {
                for (unsigned int j = 0; j < num_isvs; j++) {
                    (*dPK2StressdPreviousISVs.value)[dimension * num_isvs * i + num_isvs * i + j] +=
                        (*get_dPK2MeanStressdPreviousISVs())[j];
                }
            }
        }

        void residual::setdCauchyStressdT() {
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the temperature
             */

            auto dCauchyStressdT = get_SetDataStorage_dCauchyStressdT();

            dCauchyStressdT.zero(sot_dimension);

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dCauchyStressdPK2Stress(get_dCauchyStressdPK2Stress()->data(), sot_dimension, sot_dimension);

            Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_dPK2StressdT(get_dPK2StressdT()->data(),
                                                                                        sot_dimension);

            Eigen::Map<Eigen::Vector<floatType, sot_dimension> > map_dCauchyStressdT(dCauchyStressdT.value->data(),
                                                                                     sot_dimension);

            map_dCauchyStressdT = (map_dCauchyStressdPK2Stress * map_dPK2StressdT).eval();
        }

        void residual::setdCauchyStressdPreviousT() {
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the previous temperature
             */

            auto dCauchyStressdPreviousT = get_SetDataStorage_dCauchyStressdPreviousT();

            dCauchyStressdPreviousT.zero(sot_dimension);

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dCauchyStressdPK2Stress(get_dCauchyStressdPK2Stress()->data(), sot_dimension, sot_dimension);

            Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_dPK2StressdPreviousT(
                get_dPK2StressdPreviousT()->data(), sot_dimension);

            Eigen::Map<Eigen::Vector<floatType, sot_dimension> > map_dCauchyStressdPreviousT(
                dCauchyStressdPreviousT.value->data(), sot_dimension);

            map_dCauchyStressdPreviousT = (map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousT).eval();
        }

        void residual::setdCauchyStressdPreviousISVs() {
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the previous internal state variables
             */

            const unsigned int num_isvs = get_dPK2StressdPreviousISVs()->size() / sot_dimension;

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dCauchyStressdPK2Stress(get_dCauchyStressdPK2Stress()->data(), sot_dimension, sot_dimension);

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, -1, Eigen::RowMajor> > map_dPK2StressdPreviousISVs(
                get_dPK2StressdPreviousISVs()->data(), sot_dimension, num_isvs);

            auto dCauchyStressdPreviousISVs = get_SetDataStorage_dCauchyStressdPreviousISVs();

            dCauchyStressdPreviousISVs.zero(sot_dimension * num_isvs);

            Eigen::Map<Eigen::Matrix<floatType, sot_dimension, -1, Eigen::RowMajor> > map_dCauchyStressdPreviousISVs(
                dCauchyStressdPreviousISVs.value->data(), sot_dimension, num_isvs);

            map_dCauchyStressdPreviousISVs = (map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousISVs).eval();
        }

        void residual::setPreviousdCauchyStressdT() {
            /*!
             * Set previous the derivative of the Cauchy stress w.r.t. the temperature
             */

            auto previousdCauchyStressdT = get_SetDataStorage_previousdCauchyStressdT();

            Eigen::Map<const Eigen::Matrix<floatType, sot_dimension, sot_dimension, Eigen::RowMajor> >
                map_dCauchyStressdPK2Stress(get_dCauchyStressdPK2Stress()->data(), sot_dimension, sot_dimension);

            Eigen::Map<const Eigen::Vector<floatType, sot_dimension> > map_previousdPK2StressdT(
                get_previousdPK2StressdT()->data(), sot_dimension);

            previousdCauchyStressdT.zero(sot_dimension);

            Eigen::Map<Eigen::Vector<floatType, sot_dimension> > map_previousdCauchyStressdT(
                previousdCauchyStressdT.value->data(), sot_dimension);

            map_previousdCauchyStressdT = (map_dCauchyStressdPK2Stress * map_previousdPK2StressdT).eval();
        }

        void residual::setdRdT() {
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            auto dRdT = get_SetDataStorage_dRdT();

            *dRdT.value = *get_dCauchyStressdT();
        }

        void residual::setdRdT(const floatVector &dRdT) {
            /*!
             * Pass-through function to ResidualBase::setdRdT
             *
             * Required because of overloading
             *
             * \param &dRdT: The derivative of the residual w.r.t. the temperature
             */

            tardigradeHydra::ResidualBase<>::setdRdT(dRdT);
        }

    }  // namespace linearViscoelasticity

}  // namespace tardigradeHydra
