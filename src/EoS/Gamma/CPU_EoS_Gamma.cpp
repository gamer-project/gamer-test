#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO )



/********************************************************
1. Ideal gas EoS with a constant adiabatic index (EOS_GAMMA)

2. This file is shared by both CPU and GPU

   GPU_EoS_Gamma.cu -> CPU_EoS_Gamma.cpp

3. Three steps are required to implement an EoS

   I.   Set EoS auxiliary arrays
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions
********************************************************/



// =============================================
// I. Set EoS auxiliary arrays
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_Gamma
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
//                   AuxArray_Flt[0] = gamma
//                   AuxArray_Flt[1] = gamma-1
//                   AuxArray_Flt[2] = 1/(gamma-1)
//                   AuxArray_Flt[3] = 1/gamma
//                   AuxArray_Flt[4] = (mean molecular weight)*(atomic mass unit)/(Boltzmann constant)*(UNIT_E/UNIT_M)
//                   AuxArray_Flt[5] = 1/AuxArray_Flt[4]
//
// Note        :  1. Invoked by EoS_Init_Gamma()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Physical constants such as Const_amu/Const_kB should be set to unity when disabling OPT__UNIT
//                5. Do not change the order of AuxArray_Flt[]
//                   --> For example, the dual-energy routines assume AuxArray_Flt[0]=GAMMA
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_Gamma( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[0] = GAMMA;
   AuxArray_Flt[1] = GAMMA - 1.0;
   AuxArray_Flt[2] = 1.0 / ( GAMMA - 1.0 );
   AuxArray_Flt[3] = 1.0 / GAMMA;
   AuxArray_Flt[4] = ( OPT__UNIT ) ? MOLECULAR_WEIGHT * Const_amu / Const_kB * (UNIT_E/UNIT_M)
                                   : MOLECULAR_WEIGHT;
   AuxArray_Flt[5] = 1.0 / AuxArray_Flt[4];

} // FUNCTION : EoS_SetAuxArray_Gamma
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
//     (4) EoS_DensEint2Temp_* [OPTIONAL]
//     (5) EoS_DensTemp2Pres_* [OPTIONAL]
//     (6) EoS_General_*       [OPTIONAL]
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_Gamma
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Gamma() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//                ExtraInOut : Useless for this EoS
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_Gamma( const real Dens, const real Eint, const real Passive[],
                                     const double AuxArray_Flt[], const int AuxArray_Int[],
                                     const real *const Table[EOS_NTABLE_MAX], real ExtraInOut[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "input density",        __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Eint, NULL, "input internal energy",__FILE__, __FUNCTION__, __LINE__, true );
#  endif // GAMER_DEBUG


   const real Gamma_m1 = (real)AuxArray_Flt[1];
   real Pres;

   Pres = Eint * Gamma_m1;

   return Pres;

} // FUNCTION : EoS_DensEint2Pres_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_Gamma
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_Gamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//                ExtraInOut : Useless for this EoS
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_Gamma( const real Dens, const real Pres, const real Passive[],
                                     const double AuxArray_Flt[], const int AuxArray_Int[],
                                     const real *const Table[EOS_NTABLE_MAX], real ExtraInOut[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "input density",  __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Pres, NULL, "input pressure", __FILE__, __FUNCTION__, __LINE__, true );
#  endif // GAMER_DEBUG


   const real _Gamma_m1 = (real)AuxArray_Flt[2];
   real Eint;

   Eint = Pres * _Gamma_m1;

   return Eint;

} // FUNCTION : EoS_DensPres2Eint_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_Gamma
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_Gamma()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//                ExtraInOut : Useless for this EoS
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_Gamma( const real Dens, const real Pres, const real Passive[],
                                     const double AuxArray_Flt[], const int AuxArray_Int[],
                                     const real *const Table[EOS_NTABLE_MAX], real ExtraInOut[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "input density",  __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Pres, NULL, "input pressure", __FILE__, __FUNCTION__, __LINE__, true );
#  endif // GAMER_DEBUG


   const real Gamma = (real)AuxArray_Flt[0];
   real Cs2;

   Cs2 = Gamma * Pres / Dens;

   return Cs2;

} // FUNCTION : EoS_DensPres2CSqr_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Temp_Gamma
// Description :  Convert gas mass density and internal energy density to gas temperature
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Gamma() for the values stored in AuxArray_Flt/Int[]
//                3. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//                ExtraInOut : Useless for this EoS
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Temp_Gamma( const real Dens, const real Eint, const real Passive[],
                                     const double AuxArray_Flt[], const int AuxArray_Int[],
                                     const real *const Table[EOS_NTABLE_MAX], real ExtraInOut[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "input density",         __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Eint, NULL, "input internal energy", __FILE__, __FUNCTION__, __LINE__, true );
#  endif // GAMER_DEBUG


   const real Gamma_m1 = (real)AuxArray_Flt[1];
   const real m_kB     = (real)AuxArray_Flt[4];
   real Pres, Temp;

   Pres = Eint * Gamma_m1;
   Temp = m_kB * Pres / Dens;

   return Temp;

} // FUNCTION : EoS_DensEint2Temp_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensTemp2Pres_Gamma
// Description :  Convert gas mass density and temperature to gas pressure
//
// Note        :  1. See EoS_SetAuxArray_Gamma() for the values stored in AuxArray_Flt/Int[]
//                2. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Temp       : Gas temperature in kelvin
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//                ExtraInOut : Useless for this EoS
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensTemp2Pres_Gamma( const real Dens, const real Temp, const real Passive[],
                                     const double AuxArray_Flt[], const int AuxArray_Int[],
                                     const real *const Table[EOS_NTABLE_MAX], real ExtraInOut[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "input density",     __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Temp, NULL, "input temperature", __FILE__, __FUNCTION__, __LINE__, true );
#  endif // GAMER_DEBUG


   const real _m_kB = (real)AuxArray_Flt[5];
   real Pres;

   Pres = Temp * Dens * _m_kB;

   return Pres;

} // FUNCTION : EoS_DensTemp2Pres_Gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_General_Gamma
// Description :  General EoS converter: In[] -> Out[]
//
// Note        :  1. See EoS_DensEint2Pres_Gamma()
//                2. In[] and Out[] must NOT overlap
//                3. Useless for this EoS
//
// Parameter   :  Mode       : To support multiple modes in this general converter
//                Out        : Output array
//                In         : Input array
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_General_Gamma( const int Mode, real Out[], const real In[], const double AuxArray_Flt[],
                               const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// not used by this EoS

} // FUNCTION : EoS_General_Gamma



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Gamma;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_Gamma;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_Gamma;
FUNC_SPACE EoS_DE2T_t EoS_DensEint2Temp_Ptr = EoS_DensEint2Temp_Gamma;
FUNC_SPACE EoS_DT2P_t EoS_DensTemp2Pres_Ptr = EoS_DensTemp2Pres_Gamma;
FUNC_SPACE EoS_GENE_t EoS_General_Ptr       = EoS_General_Gamma;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_Gamma
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_Gamma()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_Gamma( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//                EoS_DensEint2Temp_CPU/GPUPtr : ...
//                EoS_DensTemp2Pres_CPU/GPUPtr : ...
//                EoS_General_CPU/GPUPtr       : ...
//
// Return      :  EoS_DensEint2Pres_CPU/GPUPtr, EoS_DensPres2Eint_CPU/GPUPtr,
//                EoS_DensPres2CSqr_CPU/GPUPtr, EoS_DensEint2Temp_CPU/GPUPtr,
//                EoS_DensTemp2Pres_CPU/GPUPtr, EoS_General_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_Gamma( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                           EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                           EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr,
                           EoS_DE2T_t &EoS_DensEint2Temp_GPUPtr,
                           EoS_DT2P_t &EoS_DensTemp2Pres_GPUPtr,
                           EoS_GENE_t &EoS_General_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Temp_GPUPtr, EoS_DensEint2Temp_Ptr, sizeof(EoS_DE2T_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensTemp2Pres_GPUPtr, EoS_DensTemp2Pres_Ptr, sizeof(EoS_DT2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_General_GPUPtr,       EoS_General_Ptr,       sizeof(EoS_GENE_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_Gamma( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                           EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                           EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr,
                           EoS_DE2T_t &EoS_DensEint2Temp_CPUPtr,
                           EoS_DT2P_t &EoS_DensTemp2Pres_CPUPtr,
                           EoS_GENE_t &EoS_General_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
   EoS_DensEint2Temp_CPUPtr = EoS_DensEint2Temp_Ptr;
   EoS_DensTemp2Pres_CPUPtr = EoS_DensTemp2Pres_Ptr;
   EoS_General_CPUPtr       = EoS_General_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_Gamma( double [], int [] );
void EoS_SetCPUFunc_Gamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_GENE_t & );
#ifdef GPU
void EoS_SetGPUFunc_Gamma( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_GENE_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_Gamma
// Description :  Initialize EoS
//
// Note        :  1. Set auxiliary arrays by invoking EoS_SetAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU EoS routines by invoking EoS_SetCPU/GPUFunc_*()
//                3. Invoked by EoS_Init()
//                   --> Enable it by linking to the function pointer "EoS_Init_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void EoS_Init_Gamma()
{

   EoS_SetAuxArray_Gamma( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_Gamma( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                         EoS_DensPres2CSqr_CPUPtr, EoS_DensEint2Temp_CPUPtr,
                         EoS_DensTemp2Pres_CPUPtr, EoS_General_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_Gamma( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr,
                         EoS_DensPres2CSqr_GPUPtr, EoS_DensEint2Temp_GPUPtr,
                         EoS_DensTemp2Pres_GPUPtr, EoS_General_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_Gamma

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO )
