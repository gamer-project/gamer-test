#ifndef __CUFLU_DUALENERGY__
#define __CUFLU_DUALENERGY__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  !defined SRHD && defined DUAL_ENERGY )



// internal functions
#if ( DUAL_ENERGY == DE_ENPY  &&  defined __CUDACC__ )
GPU_DEVICE
static real Hydro_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 );
GPU_DEVICE
static real Hydro_DensEntropy2Pres( const real Dens, const real Enpy, const real Gamma_m1,
                                    const bool CheckMinPres, const real MinPres );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DualEnergyFix
// Description :  Correct the internal and total energies using the dual-energy formalism
//
// Note        :  1. Invoked by Hydro_FullStepUpdate(), InterpolateGhostZon(), ...
//                2. A floor value "MinPres" is applied to the corrected pressure if CheckMinPres is on
//                3  A floor value "TINY_NUMBER" is applied to the input entropy as well
//                4. Call-by-reference for "Etot, Enpy, and DE_Status"
//                5. Fluid variables returned by this function are guaranteed to be consistent with each other
//                   --> They must satisfy "entropy = pressure / density^(Gamma-1)", where pressure is calculated
//                       by (Etot - Ekin - Emag)*(Gamma-1.0)
//                   --> It doesn't matter we use entropy to correct Eint or vice versa, and it also holds even when
//                       the floor value is applied to pressure
//
// Parameter   :  Dens             : Mass density
//                MomX/Y/Z         : Momentum density
//                Etot             : Total energy density
//                Enpy             : Entropy
//                DE_Status        : Assigned to (DE_UPDATED_BY_ETOT / DE_UPDATED_BY_DUAL / DE_UPDATED_BY_MIN_PRES)
//                                   to indicate whether this cell is updated by the total energy, dual energy variable,
//                                   or pressure floor (MinPres)
//                Gamma_m1         : Adiabatic index - 1.0
//                _Gamma_m1        : 1.0/Gamma_m1
//                CheckMinPres     : Return Hydro_CheckMinPres()
//                                   --> In some cases we actually want to check if pressure becomes unphysical,
//                                       for which we don't want to enable this option
//                MinPres          : Minimum allowed pressure
//                DualEnergySwitch : if ( Eint/(Ekin+Emag) < DualEnergySwitch ) ==> correct Eint and Etot
//                                   else                                       ==> correct Enpy
//                Emag             : Magnetic energy density (0.5*B^2) --> for MHD only
//
// Return      :  Etot, Enpy, DE_Status
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Enpy, char &DE_Status, const real Gamma_m1, const real _Gamma_m1,
                          const bool CheckMinPres, const real MinPres, const real DualEnergySwitch,
                          const real Emag )
{

   const bool CheckMinPres_No = false;
   const bool CheckMinEint_No = false;

// apply entropy floor
   Enpy = FMAX( Enpy, TINY_NUMBER );


// calculate energies
// --> note that here Eint can even be negative due to numerical errors
// --> Enth (i.e., non-thermal energy) includes both kinetic and magnetic energies
   real Enth, Eint, Pres;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Etot, NULL, NULL,
                          CheckMinEint_No, NULL_REAL, Emag );
   Enth = Etot - Eint;


// determine whether or not to use the dual-energy variable (entropy or internal energy) to correct the total energy density
   if ( Eint/Enth < DualEnergySwitch )
   {
//    correct total energy
//    --> we will apply pressure floor later
#     if   ( DUAL_ENERGY == DE_ENPY )
      Pres = Hydro_DensEntropy2Pres( Dens, Enpy, Gamma_m1, CheckMinPres_No, NULL_REAL );
      Eint = Pres*_Gamma_m1;

#     elif ( DUAL_ENERGY == DE_EINT )
#     error : DE_EINT is NOT supported yet !!
#     endif

      Etot      = Enth + Eint;
      DE_Status = DE_UPDATED_BY_DUAL;
   }

   else
   {
//    correct entropy
      Pres      = Eint*Gamma_m1;
      Enpy      = Hydro_DensPres2Entropy( Dens, Pres, Gamma_m1 );
      DE_Status = DE_UPDATED_BY_ETOT;
   } // if ( Eint/Enth < DualEnergySwitch ) ... else ...


// apply pressure floor
   if ( CheckMinPres  &&  Pres < MinPres )
   {
      Pres = MinPres;
      Eint = Pres*_Gamma_m1;

//    ensure that both energy and entropy are consistent with the pressure floor
      Etot      = Enth + Eint;
      Enpy      = Hydro_DensPres2Entropy( Dens, Pres, Gamma_m1 );
      DE_Status = DE_UPDATED_BY_MIN_PRES;
   }

} // FUNCTION : Hydro_DualEnergyFix



#if ( DUAL_ENERGY == DE_ENPY )

// Hydro_Con2Entropy() is used by CPU only
#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Entropy
// Description :  Evaluate the gas entropy from the input fluid variables
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_Init_ByFunction_AssignData(), Gra_Close(), Init_ByFile(), ...
//                3. Currently this function does NOT apply pressure floor when calling Hydro_Con2Pres()
//                   --> However, note that Hydro_DensPres2Entropy() does apply a floor value (TINY_NUMBER) for entropy
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Total energy density
//                Emag              : Magnetic energy density (0.5*B^2) --> for MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables
//
// Return      :  Enpy
//-------------------------------------------------------------------------------------------------------
real Hydro_Con2Entropy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                        const real Emag, const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[],
                        const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// currently this function does NOT apply pressure floor when calling Hydro_Con2Pres()
   const bool CheckMinPres_No = false;

   real Pres, Enpy;

// calculate pressure and convert it to entropy
// --> note that DE_ENPY only works with EOS_GAMMA, which does not involve passive scalars
   Pres = Hydro_Con2Pres( Dens, MomX, MomY, MomZ, Engy, NULL, CheckMinPres_No, NULL_REAL, Emag,
                          EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   Enpy = Hydro_DensPres2Entropy( Dens, Pres, EoS_AuxArray_Flt[1] );

   return Enpy;

} // FUNCTION : Hydro_Con2Entropy
#endif // ifndef __CUDACC__



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensPres2Entropy
// Description :  Evaluate the gas entropy from the input density and pressure
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_Con2Entropy() and Hydro_DualEnergyFix()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value (TINY_NUMBER) is applied to the returned value
//
// Parameter   :  Dens     : Mass density
//                Pres     : Pressure
//                Gamma_m1 : Adiabatic index - 1.0
//
// Return      :  Enpy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 )
{

   real Enpy;

// calculate entropy
   Enpy = Pres*POW( Dens, -Gamma_m1 );

// apply a floor value
   Enpy = FMAX( Enpy, TINY_NUMBER );

   return Enpy;

} // FUNCTION : Hydro_DensPres2Entropy



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensEntropy2Pres
// Description :  Evaluate the gas pressure from the input density and entropy
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_DualEnergyFix(), Flu_Close(), Hydro_Aux_Check_Negative(), and Flu_FixUp()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value "MinPres" is applied to the returned pressure if CheckMinPres is on
//
// Parameter   :  Dens         : Mass density
//                Enpy         : Enpy
//                Gamma_m1     : Adiabatic index - 1.0
//                CheckMinPres : Return Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//
// Return      :  Pres
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_DensEntropy2Pres( const real Dens, const real Enpy, const real Gamma_m1,
                             const bool CheckMinPres, const real MinPres )
{

   real Pres;

// calculate pressure
   Pres = Enpy*POW( Dens, Gamma_m1 );

// apply a floor value
   if ( CheckMinPres )  Pres = Hydro_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : Hydro_DensEntropy2Pres

#endif // #if ( DUAL_ENERGY == DE_ENPY )



#endif // #if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )



#endif // #ifndef __CUFLU_DUALENERGY__
