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

    //! The base class for deformations
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

            secondOrderTensor getSubConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex);

            secondOrderTensor getPrecedingConfiguration(const unsigned int &index);

            secondOrderTensor getFollowingConfiguration(const unsigned int &index);

            secondOrderTensor getConfiguration(const unsigned int &index);

            secondOrderTensor getPreviousSubConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex);

            secondOrderTensor getPreviousPrecedingConfiguration(const unsigned int &index);

            secondOrderTensor getPreviousFollowingConfiguration(const unsigned int &index);

            secondOrderTensor getPreviousConfiguration(const unsigned int &index);

            floatVector getSubConfigurationJacobian(const unsigned int &lowerIndex, const unsigned int &upperIndex);

            floatVector getPrecedingConfigurationJacobian(const unsigned int &index);

            floatVector getFollowingConfigurationJacobian(const unsigned int &index);

            // CACHED DATA STORAGE OPERATIONS
            virtual void addIterationData(dataBase *data) override;

            virtual void addNLStepData(dataBase *data) override;
            // END CACHED DATA STORAGE OPERATIONS

            // PASS-THROUGH FUNCTIONS
            unsigned int getNumConfigurations(); //TODO: Extract this from hydra to here
            // END PASS-THROUGH FUNCTIONS

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

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, inverseConfigurations, floatVector, passThrough)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousInverseConfigurations, floatVector, passThrough)

    };

}

#include "tardigrade_DeformationBase.tpp"

#endif
