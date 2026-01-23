/**
 ******************************************************************************
 * \file tardigrade_CoreDefinitions.h
 ******************************************************************************
 * Core definitions for tardigrade hydra
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_COREDEFINITIONS
#define TARDIGRADE_HYDRA_COREDEFINITIONS

#include <vector>

namespace tardigradeHydra {

    // forward class declarations
    class hydraBase;  //!< The base hydra class
    template <class container = hydraBase>
    class ResidualBase;        //!< The base residual class
    class SolverBase;          //!< The base solver class
    class SolverStepBase;      //!< The base solver step class
    class TrialStepBase;       //!< The base trial-step class
    class PreconditionerBase;  //!< The base preconditioner class

    namespace unit_test {
        class hydraBaseTester;           //!< Friend class for hydraBase for unit testing
        class SolverStepBaseTester;      //!< Friend class for SolverStepBase for unit testing
        class SolverBaseTester;          //!< Friend class for SolverBase for unit testing
        class PreconditionerBaseTester;  //!< Friend class for PreconditionerBase for unit testing
    }  // namespace unit_test

    typedef double                               floatType;    //!< Define the float values type.
    typedef std::vector<floatType>               floatVector;  //!< Define a vector of floats
    typedef std::vector<std::vector<floatType> > floatMatrix;  //!< Define a matrix of floats

    // Define tensors of known size
    typedef std::vector<floatType> dimVector;          //!< Dimension vector
    typedef std::vector<floatType> secondOrderTensor;  //!< Second order tensors
    typedef std::vector<floatType> thirdOrderTensor;   //!< Third order tensors
    typedef std::vector<floatType> fourthOrderTensor;  //!< Fourth order tensors

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
    typedef libxsmm_mmfunction<floatType> kernel_type;  //!< The libxsmm kernel type
#endif

    typedef void (hydraBase::*hydraBaseFxn)();  //!< Typedef for passing pointers to hydraBase functions

}  // namespace tardigradeHydra

#include "tardigrade_CoreDefinitions.cpp"

#endif
