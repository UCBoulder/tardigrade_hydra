/**
  ******************************************************************************
  * \file tardigrade_CoreDefinitions.h
  ******************************************************************************
  * Core definitions for tardigrade hydra
  ******************************************************************************
  */ 

#ifndef TARDIGRADE_HYDRA_COREDEFINITIONS
#define TARDIGRADE_HYDRA_COREDEFINITIONS

namespace tardigradeHydra{

    // forward class declarations
    class hydraBase; //!< The base hydra class
    class residualBase; //!< The base residual class

    namespace unit_test{
        class hydraBaseTester; //!< Friend class for hydraBase for unit testing
    }

}

#include"tardigrade_CoreDefinitions.cpp"

#endif
