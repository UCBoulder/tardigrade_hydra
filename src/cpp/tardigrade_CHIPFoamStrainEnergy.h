/**
 ******************************************************************************
 * \file tardigrade_CHIPFoamStrainEnergy.h
 ******************************************************************************
 * The CHIPFoam strain-energy potential
 ******************************************************************************
 */

#ifndef TARDIGRADE_CHIPFOAMSTRAINENERGY
#define TARDIGRADE_CHIPFOAMSTRAINENERGY

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

    /*!
     * The CHIPFoam strain-energy potential
     */
    class CHIPFoamStrainEnergy : public HyperelasticBase {
       public:
        /*!
         * Default residual
         */
        CHIPFoamStrainEnergy() : HyperelasticBase(), _parameters({}) {};

        /*!
         * Main utilization constructor
         *
         * \param *_hydra: A pointer to a hydraBase object
         * \param &_numEquations: The number of equations defined by the residual
         * \param &parameters: The parameter vector organized as
         *    Khat, Ghat, Jb
         */
        CHIPFoamStrainEnergy(hydraBase *_hydra, const unsigned int &_numEquations, const floatVector &parameters)
            : HyperelasticBase(_hydra, _numEquations), _parameters(parameters){

            TARDIGRADE_ERROR_TOOLS_CHECK(_parameters.size() == 3, "The parameters vector must have a size of 3");

            setInitialized();

        }

        const floatType get_Khat();

        const floatType get_Ghat();

        const floatType get_Jb();

       protected:
        //! The model parameters
        floatVector _parameters;

        virtual void setJe(bool isPrevious);

        virtual void setJe();

        virtual void setPreviousJe();

        virtual void setdJedFe(bool isPrevious);

        virtual void setdJedFe();

        virtual void setdPreviousJedPreviousFe();

        virtual void setd2JedFe2(bool isPrevious);

        virtual void setd2JedFe2();

        virtual void setd2PreviousJedPreviousFe2();

        virtual void setIbar1(bool isPrevious);

        virtual void setIbar1();

        virtual void setPreviousIbar1();

        virtual void setdIbar1dFe(bool isPrevious);

        virtual void setdIbar1dFe();

        virtual void setdPreviousIbar1dPreviousFe();

        virtual void setd2Ibar1dFe2(bool isPrevious);

        virtual void setd2Ibar1dFe2();

        virtual void setd2PreviousIbar1dPreviousFe2();

        virtual void setWLB(bool isPrevious);

        virtual void setWLB();

        virtual void setPreviousWLB();

        //! Check if the class has been initialized
        const bool isInitialized(){ return is_initialized; }

        //! Set that the class has been initialized
        void setInitialized(){ is_initialized = true; };

       private:
        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Je, floatType, setJe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousJe, floatType, setPreviousJe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dJedFe, secondOrderTensor, setdJedFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousJedPreviousFe, secondOrderTensor, setdPreviousJedPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2JedFe2, fourthOrderTensor, setd2JedFe2)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousJedPreviousFe2, fourthOrderTensor, setd2PreviousJedPreviousFe2)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Ibar1, floatType, setIbar1)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousIbar1, floatType, setPreviousIbar1)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dIbar1dFe, secondOrderTensor, setdIbar1dFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousIbar1dPreviousFe, secondOrderTensor, setdPreviousIbar1dPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2Ibar1dFe2, fourthOrderTensor, setd2Ibar1dFe2)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousIbar1dPreviousFe2, fourthOrderTensor, setd2PreviousIbar1dPreviousFe2)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, WLB, floatType, setWLB);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousWLB, floatType, setPreviousWLB);

        //! Whether the class has been initialized or not
        bool is_initialized = false;

    };

}  // namespace tardigradeHydra

#include "tardigrade_CHIPFoamStrainEnergy.tpp"

#endif
