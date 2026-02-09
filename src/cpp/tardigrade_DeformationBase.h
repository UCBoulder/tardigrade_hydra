/**
 ******************************************************************************
 * \file tardigrade_DeformationBase.h
 ******************************************************************************
 * The base class for defining the material deformation
 ******************************************************************************
 */

#ifndef TARDIGRADE_DEFORMATIONBASE_H
#define TARDIGRADE_DEFORMATIONBASE_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_SetDataStorage.h"

namespace tardigradeHydra {

    namespace unit_test{

        class DeformationBaseTester; //!< Frient class for DeformationBase unit testing

    }

    class DeformationBase : public CachingDataBase {

        public:

            /*!
             * Constructor of DeformationBase
             *
             * \param *_hydra: The containing hydra object
             */
            DeformationBase(hydraBase *_hydra=nullptr) : hydra(_hydra) { }

            template<unsigned int leading_rows, unsigned int size, unsigned int dim>
            floatVector getSubConfiguration(const floatVector &configurations, const unsigned int &lowerIndex,
                                            const unsigned int &upperIndex);

            template<unsigned int leading_rows, unsigned int size, unsigned int dim>
            floatVector getSubConfigurationJacobian(const floatVector  &configurations,
                                                    const unsigned int &lowerIndex, const unsigned int &upperIndex);

            // CACHED DATA STORAGE OPERATIONS
            virtual void addIterationData(dataBase *data) override;

            virtual void addNLStepData(dataBase *data) override;
            // END CACHED DATA STORAGE OPERATIONS

            //! Pointer to the containing hydraBase object
            hydraBase *hydra;

        protected:

        private:
            friend class unit_test::DeformationBaseTester;  //!< Friend class which allows modification of private variables. ONLY
                                                            //!< TO BE USED FOR TESTING!

            //! A pass through function that does nothing
            void passThrough() {}

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, configurations, floatVector, passThrough)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousConfigurations, floatVector, passThrough)

    };

}

#include "tardigrade_DeformationBase.tpp"

#endif
