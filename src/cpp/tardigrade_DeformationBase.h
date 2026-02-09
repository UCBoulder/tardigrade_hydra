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

            template<unsigned int leading_rows, unsigned int size, unsigned int dim>
            floatVector getSubConfiguration(const floatVector &configurations, const unsigned int &lowerIndex,
                                            const unsigned int &upperIndex);

        protected:

        private:

    };

}

#include "tardigrade_DeformationBase.cpp"

#endif
