/**
 ******************************************************************************
 * \file tardigrade_NeoHookianStrainEnergy.h
 ******************************************************************************
 * A Neo-Hookian strain energy function
 ******************************************************************************
 */

#ifndef TARDIGRADE_NEOHOOKIANSTRAINENERGY
#define TARDIGRADE_NEOHOOKIANSTRAINENERGY

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

    /*!
     * A Neo-Hookian strain energy
     */
    class NeoHookianStrainEnergy : public HyperelasticBase {
       public:
        /*!
         * Default constructor
         */
        NeoHookianStrainEnergy() : HyperelasticBase(), _C10(0), _D1(0) {};

        /*!
         * Main utilization constructor
         *
         * \param *_hydra: A pointer to a hydraBase object
         * \param &_numEquations: The number of equations defined by the residual
         * \param &parameters: The parameter vector organized as
         *    C10, D1
         */
        NeoHookianStrainEnergy(hydraBase *_hydra, const unsigned int &_numEquations, const floatVector &parameters)
            : HyperelasticBase(_hydra, _numEquations) {
            TARDIGRADE_ERROR_TOOLS_CHECK(parameters.size() == 2, "The parameters vector must have a size of 2");

            _C10 = parameters[0];
            _D1  = parameters[1];

            setInitialized();
        }

        virtual void setStrainEnergy(const bool isPrevious) override;

        virtual void setStrainEnergyJacobians(const bool isPrevious) override;

        virtual void setStrainEnergyHessians(const bool isPrevious) override;

       protected:
        //! Check if the class has been initialized
        const bool isInitialized() { return is_initialized; }

        //! Set that the class has been initialized
        void setInitialized() { is_initialized = true; };

        //! The parameter associated with Ibar1
        floatType _C10;

        //! The parameter associated with the volumetric deformation
        floatType _D1;

       private:
        //! Whether the class has been initialized or not
        bool is_initialized = false;
    };

}  // namespace tardigradeHydra

#include "tardigrade_NeoHookianStrainEnergy.tpp"

#endif
