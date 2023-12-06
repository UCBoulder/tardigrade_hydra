/**
  ******************************************************************************
  * \file tardigrade_hydra_umat.cpp
  ******************************************************************************
  * The Abaqus UMAT c++ interface template
  ******************************************************************************
  */

#include<tardigrade_hydra_umat.h>

extern "C" void umat_( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                       double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                       const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                       const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                       const int &NSHR,      const int &NTENS,     const int &NSTATEV,    const double *PROPS,  const int &NPROPS,
                       const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                       const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                       const int *JSTEP,     const int &KINC ){
    /*!
     * A template Abaqus UMAT c++ interface.
     *
     * The variables defined in this interface are described more completely in the Abaqus User Subroutines
     * manual entry for UMAT user subroutines.
     *
     * The 80 character assumption for string conversion is based on the data line limits on character strings as
     * described in the Abaqus User's Manual > Introduction & Spatial Modeling > Input Syntax Rules and the sample
     * Fortran UMAT in the Abaqus User's Manual > User Subroutines > Abaqus/Standard User Subroutines > UMAT.
     *
     * \param *STRESS: Cauchy stress tensor at beginning of time increment stored in vector form, \f$ \sigma \f$.
     * \param *STATEV: State variable vector.
     * \param *DDSDDE: Jacobian matrix \f$ \frac{\delta \Delta \sigma}{\delta \Delta \epsilon} \f$.
     * \param &SSE: Specific elastic strain energy.
     * \param &SPD: Specific plastic dissipation.
     * \param &SCD: Specific "creep" dissipation.
     * \param &RPL: Volumetric heat generation per unit time at the end of the time increment.
     * \param *DDSDDT: Variation of stress increment with respect to temperature.
     * \param *DRPLDE: Variation of RPL with respect to strain increment.
     * \param &DRPLDT: Variation of RPL with respect to temperature .
     * \param *STRAN: Strain tensor at the beginning of the time increment stored in vector form, \f$ \epsilon \f$.
     * \param *DSTRAN: Strain increment tensor stored in vector form.
     * \param *TIME: Time vector at beginning of time increment = {Step time, Total time}.
     * \param &DTIME: Time increment.
     * \param &TEMP: Temperature at the start of the time increment.
     * \param &DTEMP: Temperature increment.
     * \param *PREDEF: Interpolated predefined field variables at current integration point.
     * \param *DPRED: Change in predefined field variables.
     * \param *CMNAME: User defined material name. Left justified as passed by FORTRAN.
     * \param &NDI: Number of direct stress components at this integration point.
     * \param &NSHR: Number of engineering shear stress components at this integration point.
     * \param &NTENS: Size of the stress and strain component array. NTENS = NDI + NSHR
     * \param &NSTATEV: Number of state variables for this material, CMNAME.
     * \param *PROPS: Material model constants defined as part of the *MATERIAL keyword in the input file.
     * \param &NPROPS: Number of user defined material constants.
     * \param *COORDS: Coordinates of the current Gauss point.
     * \param *DROT: Rigid body rotation increment matrix.
     * \param &PNEWDT: Ratio of suggested new time increment to current time increment, e.g. the next increment's DTIME.
     * \param &CELENT: Characteristic element length.
     * \param *DFGRD0: Deformation gradient matrix at the beginning of the time increment.
     * \param *DFGRD1: Deformation gradient at the end of the current time increment.
     * \param &NOEL: Current element number.
     * \param &NPT: Integration point number.
     * \param &LAYER: Layer number for composite shells and layered solids.
     * \param &KSPT: Section point number within the current layer.
     * \param *JSTEP: Meta data vector = {Step number, procedure type key, NLGEOM switch, linear perturbation switch}
     * \param &KINC: Increment number.
     */

     //Add switching logic to handle more than one UMAT.
     //Call the appropriate UMAT interface
     tardigradeHydra::abaqusInterface( STRESS, STATEV, DDSDDE,    SSE,    SPD,
                                           SCD,    RPL, DDSDDT, DRPLDE, DRPLDT,
                                         STRAN, DSTRAN,   TIME,  DTIME,   TEMP,
                                         DTEMP, PREDEF,  DPRED, CMNAME,    NDI,
                                          NSHR,  NTENS, NSTATEV,  PROPS, NPROPS,
                                        COORDS,   DROT, PNEWDT, CELENT, DFGRD0,
                                        DFGRD1,   NOEL,    NPT,  LAYER,   KSPT,
                                         JSTEP,   KINC );

     return;
}
