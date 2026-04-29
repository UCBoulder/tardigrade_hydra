/**
 ******************************************************************************
 * \file tardigrade_MooneyRivlinStrainEnergy.h
 ******************************************************************************
 * A Mooney-Rivlin strain energy function
 ******************************************************************************
 */

#ifndef TARDIGRADE_MOONYRIVLINSTRAINENERGY
#define TARDIGRADE_MOONYRIVLINSTRAINENERGY

#include "tardigrade_NeoHookianStrainEnergy.h"

namespace tardigradeHydra {

    /*!
     * A Neo-Hookian strain energy
     */
    class MooneyRivlinStrainEnergy : public NeoHookianStrainEnergy {
       public:
        /*!
         * Default constructor
         */
        MooneyRivlinStrainEnergy() : NeoHookianStrainEnergy(), _C01(0) {};

        /*!
         * Main utilization constructor
         *
         * \param *_hydra: A pointer to a hydraBase object
         * \param &_numEquations: The number of equations defined by the residual
         * \param &parameters: The parameter vector organized as
         *    C10, C01, D1
         */
        MooneyRivlinStrainEnergy(hydraBase *_hydra, const unsigned int &_numEquations, const floatVector &parameters)
            : NeoHookianStrainEnergy(_hydra, _numEquations, {parameters[0], parameters[2]}) {
            TARDIGRADE_ERROR_TOOLS_CHECK(parameters.size() == 3, "The parameters vector must have a size of 3");

            _C10 = parameters[0];
            _C01 = parameters[1];
            _D1  = parameters[2];

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

        //! The parameter associated with Ibar2
        floatType _C01;

        //! The parameter associated with the volumetric deformation
        floatType _D1;

       private:
        //! Whether the class has been initialized or not
        bool is_initialized = false;
    };

}  // namespace tardigradeHydra

#include "tardigrade_MooneyRivlinStrainEnergy.tpp"

#endif
