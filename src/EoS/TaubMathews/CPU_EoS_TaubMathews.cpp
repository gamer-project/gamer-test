#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO && defined SRHD )

real VectorDotProduct( real V1, real V2, real V3 );

/********************************************************
1. Approximately relativistic ideal gas EoS (Taub-Mathews EoS)

2. This file is shared by both CPU and GPU

   GPU_EoS_TaubMathews.cu -> CPU_EoS_TaubMathews.cpp

3. Three steps are required to implement an EoS

   I.   Set an EoS auxiliary array
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions
********************************************************/



// =============================================
// I. Set an EoS auxiliary array
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_TaubMathews
// Description :  Set the auxiliary array AuxArray[]
// Note        :  1. Invoked by EoS_Init_TaubMathews()
//                2. AuxArray[] has the size of EOS_NAUX_MAX defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Do not change the order of AuxArray[]
//                   --> For example, the dual-energy routines assume AuxArray[0]=GAMMA
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_TaubMathews( double AuxArray_Flt[], int AuxArray_Int[] )
{
  return;
} // FUNCTION : EoS_SetAuxArray_TaubMathews
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_HTilde2Temp_*
//     (2) EoS_Temp2HTilde_*
//     (3) EoS_Temper2CSqr_*
// =============================================



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_GuessHTilde_TaubMathews
// Description :  Guess a reduced enthalpy for Newton-Raphson iteration
// Note        :
// Parameter   :  Con      : Conservative variables
//                         : Pointer to store guessed value
// Return      :  Constant
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_GuessHTilde_TaubMathews( const real Con[], real* const Constant, const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{
  real GuessHTilde, Discrimination;

  real Msqr = VectorDotProduct(Con[1], Con[2], Con[3]);
  real Dsqr = SQR(Con[0]);
  real abc = (real)1.0 / Dsqr;
  real E_D = Con[4] / Con[0];
  real M_Dsqr = abc * Msqr;
  real M_D = SQRT( M_Dsqr );


  // (x+y)(x-y) is more accurate than x**2-y**2
  //Constant = (E_D + M_D) * (E_D - M_D) + (real)2.0 * E_D;
  real X = SQRT( E_D*E_D + (real)2.0*E_D );
  real Y = X + M_D;
  real Z = X - M_D;
  *Constant = Y * Z;


  real A = (real)437.0 * M_Dsqr + (real)117.0;
  real B = (real)1.0 + M_Dsqr;
  real C = (real)43.0*M_Dsqr + (real)63.0;

  Discrimination  = (real)3240000.0 * SQR( B );
  Discrimination /= SQR( A );


  if ( *Constant >= Discrimination )
  {
     GuessHTilde = (real)1.3333333 * SQRT( *Constant );
  }
  else
  {
     GuessHTilde  = (real)11.18 * SQRT( B*C* (*Constant) + (real)45.0*B*B );
	 GuessHTilde -= (real)75.0 * B;
	 GuessHTilde /= C;
  }

  return GuessHTilde;
}
//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_HTilde2Temp_TaubMathews
// Description :  Convert reduced enthalpy to temperature
//
// Note        :
//
// Parameter   :  HTilde   : Reduced enthalpy
//                Temp     : Pointer to store temperature (kT/mc**)
//                DiffTemp : Pointer to store the derivitive of a reduced enthalpy with respect to temperature
//                Passive  : Passive scalars (must not used here)
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Temp, DiffTemp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_HTilde2Temp_TaubMathews( const real HTilde, real* const Temp, real* const DiffTemp,
                                         const real Passive[], const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{
  real HTildeSqr = SQR(HTilde);
  real Factor0 = (real)2.0*HTildeSqr + (real)4.0*HTilde;
  real Factor1 = SQRT( (real)9.0*HTildeSqr + (real)18.0*HTilde + (real)25.0 );
  real Factor2 = (real)5.0*HTilde + (real)5.0 + Factor1;

  if ( Temp != NULL )
    *Temp = Factor0 / Factor2;

  if ( DiffTemp != NULL )
    *DiffTemp = ( (real)4.0*HTilde + (real)4.0 ) / Factor2
                   -  Factor0 * ( ( (real)9.0*HTilde + (real)9.0 ) /  Factor1 + (real)5.0 ) / SQR(Factor2);

} // FUNCTION : EoS_HTilde2Temp_TaubMathews



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Temp2HTilde_TaubMathews
// Description :
// Note        :
// Parameter   :  Temp     : Temperature (kT/mc**)
//                Passive  : Passive scalars (must not used here)
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Reduced energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_Temp2HTilde_TaubMathews( const real Temp, const real Passive[], const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

   real HTilde; // Reduced enthalpy
   real TempSqr = Temp*Temp;
   real Factor = (real)2.25*TempSqr;

   HTilde  = Factor;
   HTilde /= (real)1.0 + SQRT( Factor + (real)1.0 );
   HTilde += (real)2.5*Temp;

   return HTilde;

} // FUNCTION : EoS_Temp2HTilde_TaubMathews



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Temper2CSqr_TaubMathews
// Description :  Convert gas temperature to sound speed squared
//
// Note        :
//
// Parameter   :  Rho      : Gas proper mass density
//                Pres     : Gas pressure
//                Passive  : Passive scalars (must not used here)
//                AuxArray : Auxiliary array (see the Note above)
//
// Return      :  Sound speed square
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_Temper2CSqr_TaubMathews( const real Rho, const real Pres, const real Passive[], const double AuxArray_Flt[],
                                         const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray == NULL )    printf( "ERROR : AuxArray == NULL in %s !!\n", __FUNCTION__ );
#  endif // GAMER_DEBUG

   real Cs2, Temp;

   Temp = Pres/Rho;

   real factor = SQRT( (real)2.25*Temp*Temp + (real)1.0 );

   Cs2  = (real) 4.5*SQR(Temp) + (real) 5.0 * Temp * factor;
   Cs2 /= (real)18.0*SQR(Temp) + (real)12.0 * Temp * factor + (real)3.0;


   return Cs2;

} // FUNCTION : EoS_Temper2CSqr_TaubMathews



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_GUESS_t  EoS_GuessHTilde_Ptr = EoS_GuessHTilde_TaubMathews;
FUNC_SPACE EoS_H2TEM_t  EoS_HTilde2Temp_Ptr = EoS_HTilde2Temp_TaubMathews;
FUNC_SPACE EoS_TEM2H_t  EoS_Temp2HTilde_Ptr = EoS_Temp2HTilde_TaubMathews;
FUNC_SPACE EoS_TEM2C_t  EoS_Temper2CSqr_Ptr = EoS_Temper2CSqr_TaubMathews;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_InitCPU/GPUFunc_TaubMathews
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_TaubMathews()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_InitFunc_TaubMathews( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_HTilde2Temp_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_Temp2HTilde_CPU/GPUPtr : ...
//                EoS_Temper2CSqr_CPU/GPUPtr : ...
//
// Return      :  EoS_HTilde2Temp_CPU, EoS_Temp2HTilde_CPU/GPUPtr, EoS_Temper2CSqr_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_TaubMathews( EoS_GUESS_t &EoS_GuessHTilde_GPUPtr,
                                 EoS_H2TEM_t &EoS_HTilde2Temp_GPUPtr,
                                 EoS_TEM2H_t &EoS_Temp2HTilde_GPUPtr,
                                 EoS_TEM2C_t &EoS_Temper2CSqr_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_GuessHTilde_GPUPtr, EoS_GuessHTilde_Ptr, sizeof(EoS_GUESS_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_HTilde2Temp_GPUPtr, EoS_HTilde2Temp_Ptr, sizeof(EoS_H2TEM_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_Temp2HTilde_GPUPtr, EoS_Temp2HTilde_Ptr, sizeof(EoS_TEM2H_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_Temper2CSqr_GPUPtr, EoS_Temper2CSqr_Ptr, sizeof(EoS_TEM2C_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_TaubMathews( EoS_GUESS_t &EoS_GuessHTilde_CPUPtr,
                                 EoS_H2TEM_t &EoS_HTilde2Temp_CPUPtr,
                                 EoS_TEM2H_t &EoS_Temp2HTilde_CPUPtr,
                                 EoS_TEM2C_t &EoS_Temper2CSqr_CPUPtr )
{
   EoS_GuessHTilde_CPUPtr = EoS_GuessHTilde_Ptr;
   EoS_HTilde2Temp_CPUPtr = EoS_HTilde2Temp_Ptr;
   EoS_Temp2HTilde_CPUPtr = EoS_Temp2HTilde_Ptr;
   EoS_Temper2CSqr_CPUPtr = EoS_Temper2CSqr_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_TaubMathews( double [] );
void EoS_SetCPUFunc_TaubMathews(EoS_GUESS_t &, EoS_H2TEM_t &, EoS_TEM2H_t &, EoS_TEM2C_t & );
#ifdef GPU
void EoS_SetGPUFunc_TaubMathews(EoS_GUESS_t &, EoS_H2TEM_t &, EoS_TEM2H_t &, EoS_TEM2C_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_TaubMathews
// Description :  Initialize EoS
//
// Note        :  1. Set an auxiliary array by invoking EoS_SetAuxArray_*()
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
void EoS_Init_TaubMathews()
{

   EoS_SetAuxArray_TaubMathews( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_TaubMathews( EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_Temper2CSqr_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_TaubMathews( EoS_GuessHTilde_GPUPtr, EoS_HTilde2Temp_GPUPtr, EoS_Temp2HTilde_GPUPtr, EoS_Temper2CSqr_GPUPtr );
#  endif

} // FUNCTION : EoS_Init_TaubMathews

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO && defined SRHD )
