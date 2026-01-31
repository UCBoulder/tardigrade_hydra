/**
 ******************************************************************************
 * \file tardigrade_PreconditionerBase.h
 ******************************************************************************
 * A C++ library for the base classes for solvers
 ******************************************************************************
 */

#ifndef TARDIGRADE_MAXROWPRECONDITIONER_H
#define TARDIGRADE_MAXROWPRECONDITIONER_H

#include "tardigrade_PreconditionerBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class MaxRowPreconditionerTester;

    }

    /*!
     * The base class for preconditioners to be used in tardigrade hydra solves
     */
    class MaxRowPreconditioner : public PreconditionerBase {
       public:
        using tardigradeHydra::PreconditionerBase::PreconditionerBase;

        virtual void preconditionVector(const floatVector &X, floatVector &Y) override;

        virtual void preconditionMatrix(const floatVector &A, floatVector &B) override;

       protected:
        void formMaxRowPreconditioner();

        virtual void formPreconditioner() override;
    };
}  // namespace tardigradeHydra

#endif
