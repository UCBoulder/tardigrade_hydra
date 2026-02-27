/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.cpp
 ******************************************************************************
 * The base class for Hyperelastic materials
 ******************************************************************************
 */

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

        /*!
         * Set the value of the elastic deformation gradient
         *
         * \param isPrevious: Flag for whether to set the current (false) or previous (true) elastic deformation
         * gradient
         */
        void HyperelasticBase::setFe(const bool isPrevious) {

            constexpr unsigned int dim     = 3;
            constexpr unsigned int sot_dim = dim * dim;

            if (isPrevious) {
                auto previousFe = get_SetDataStorage_previousFe();

                *previousFe.value =
                    secondOrderTensor(hydra->deformation->get_previousConfigurations()->begin(),
                                      hydra->deformation->get_previousConfigurations()->begin() + sot_dim);

            } else {
                auto Fe = get_SetDataStorage_Fe();

                *Fe.value = secondOrderTensor(hydra->deformation->get_configurations()->begin(),
                                              hydra->deformation->get_configurations()->begin() + sot_dim);
            }
        }

        /*!
         * Set the value of the elastic deformation gradient
         */
        void HyperelasticBase::setFe() {

            setFe(false);
        }

        /*!
         * Set the value of the previous elastic deformation gradient
         */
        void HyperelasticBase::setPreviousFe() {

            setFe(true);
        }

        /*!
         * Set the value of the derivatives of the elastic strain
         *
         * \param isPrevious: Whether this is the previous or current elastic deforamtion
         */
        void HyperelasticBase::setFeDerivatives(const bool isPrevious) {

            if (isPrevious) {
                auto previousdFedF = get_SetDataStorage_previousdFedF();

                auto previousdFedFn = get_SetDataStorage_previousdFedFn();

                *previousdFedF.value = *hydra->deformation->get_previousdF1dF();

                *previousdFedFn.value = *hydra->deformation->get_previousdF1dFn();

            } else {
                auto dFedF = get_SetDataStorage_dFedF();

                auto dFedFn = get_SetDataStorage_dFedFn();

                *dFedF.value = *hydra->deformation->get_dF1dF();

                *dFedFn.value = *hydra->deformation->get_dF1dFn();
            }
        }

        /*!
         * Set the value of the derivative of the elastic deformation gradient w.r.t. the deformation gradient
         */
        void HyperelasticBase::setdFedF() {

            setFeDerivatives(false);
        }

        /*!
         * Set the value of the derivative of the elastic deformation gradient w.r.t. the sub-deformation gradients
         */
        void HyperelasticBase::setdFedFn() {

            setFeDerivatives(false);
        }

        /*!
         * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the deformation
         * gradient
         */
        void HyperelasticBase::setPreviousdFedF() {

            setFeDerivatives(true);
        }

        /*!
         * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the sub-deformation
         * gradients
         */
        void HyperelasticBase::setPreviousdFedFn() {

            setFeDerivatives(true);
        }

        /*!
         * Set the strain energy
         *
         * \param isPrevious: A flag for whether to set the current (false) or previous (true) strain energy
         */ 
        void HyperelasticBase::setStrainEnergy(const bool isPrevious){

            throw std::runtime_error("The strain energy calculation must be defined");

        }

        /*!
         * Compute the current strain energy
         */
        void HyperelasticBase::setStrainEnergy(){

            setStrainEnergy(false);

        }

        /*!
         * Compute the previous strain energy
         */
        void HyperelasticBase::setPreviousStrainEnergy(){

            setStrainEnergy(true);

        }

        /*!
         * Set the Jacobians of the strain energy
         * 
         * \param isPrevious: Whether to set the current (false) or previous (true) Jacobians of the strain energy
         */
        void HyperelasticBase::setStrainEnergyJacobians(const bool isPrevious){

            throw std::runtime_error("The Jacobians of the strain energy must be defined");

        }

        /*!
         * Set the Jacobians of the strain energy
         */
        void HyperelasticBase::setStrainEnergyJacobians(){

            setStrainEnergyJacobians(false);

        }

        /*!
         * Set the previous Jacobians of the strain energy
         */
        void HyperelasticBase::setPreviousStrainEnergyJacobians(){

            setStrainEnergyJacobians(true);

        }

        /*!
         * Set the Hessians of the strain energy
         * 
         * \param isPrevious: Whether to set the current (false) or previous (true) Hessians of the strain energy
         */
        void HyperelasticBase::setStrainEnergyHessians(const bool isPrevious){

            throw std::runtime_error("The Hessians of the strain energy must be defined");

        }

        /*!
         * Set the Hessians of the strain energy
         */
        void HyperelasticBase::setStrainEnergyHessians(){

            setStrainEnergyHessians(false);

        }

        /*!
         * Set the previous Hessians of the strain energy
         */
        void HyperelasticBase::setPreviousStrainEnergyHessians(){

            setStrainEnergyHessians(true);

        }

        /*!
         * Set the Cauchy stress from the strain-energy function
         * 
         * \param isPrevious: Flag for whether to set the current (false) or previous (true) strain energy function
         */
        void HyperelasticBase::setCauchyStress(const bool isPrevious){

            auto dim = hydra->deformation->dimension;

            const secondOrderTensor *Fe;

            const secondOrderTensor *dStrainEnergydFe;

            SetDataStorageBase<secondOrderTensor> cauchyStress;

            if ( isPrevious ){

                Fe = get_previousFe();

                dStrainEnergydFe = get_dPreviousStrainEnergydPreviousFe();

                cauchyStress = get_SetDataStorage_previousCauchyStress();

            }
            else{

                Fe = get_Fe();

                dStrainEnergydFe = get_dStrainEnergydFe();

                cauchyStress = get_SetDataStorage_cauchyStress();

            }

            TARDIGRADE_ERROR_TOOLS_CHECK( (std::end( *Fe ) - std::begin( *Fe )) == dim * dim, "The elastic deformation must have a size of " + std::to_string( dim * dim ) );

            TARDIGRADE_ERROR_TOOLS_CHECK( (std::end( *dStrainEnergydFe ) - std::begin( *dStrainEnergydFe )) == dim * dim, "The derivative of the strain energy with respect to the elastic deformation must have a size of " + std::to_string( dim * dim ) );

            cauchyStress.zero(dim*dim);
            
            Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> Femat(Fe->data(), dim, dim); //TODO: Change this to a constant size when possible
            auto Je = Femat.determinant();

            for ( unsigned int i = 0; i < dim; ++i ){

                for ( unsigned int j = 0; j < dim; ++j ){

                    for ( unsigned int I = 0; I < dim; ++I ){

                        (*cauchyStress.value)[dim * i + j] += (*dStrainEnergydFe)[dim * i + I] * (*Fe)[dim * j + I];

                    }

                    (*cauchyStress.value)[dim * i + j] /= Je;

                }

            }

        }

        /*!
         * Set the current Cauchy stress
         */
        void HyperelasticBase::setCauchyStress(){

            setCauchyStress(false);

        }

        /*!
         * Set the previous Cauchy stress
         */
        void HyperelasticBase::setPreviousCauchyStress(){

            setCauchyStress(true);

        }

}
