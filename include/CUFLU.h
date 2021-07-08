#ifndef __CUFLU_H__
#define __CUFLU_H__



// *********************************************************************
// ** This header will be included by all CPU/GPU fluid/ELBDM solvers **
// *********************************************************************


// include "Macro.h" and "Typedef.h" here since the header "GAMER.h" is NOT included in GPU solvers
#ifdef __CUDACC__
# include "Macro.h"
# include "Typedef.h"
# include "SrcTerms.h"
# include "EoS.h"
#else
# include "GAMER.h"
#endif


// allow GPU to output messages in the debug mode
#ifdef GAMER_DEBUG
#  include "stdio.h"
#endif


// faster integer multiplication in Fermi
#if ( defined __CUDACC__  &&  __CUDA_ARCH__ >= 200 )
#  define __umul24( a, b )   ( (a)*(b) )
#  define  __mul24( a, b )   ( (a)*(b) )
#endif



// #################################
// ## macros for different models ##
// #################################

// 1. hydro macro
//=========================================================================================
#if   ( MODEL == HYDRO )

// size of different arrays
// ** to reduce the GPU memory consumption, large arrays in the fluid solvers are reused as much as possible
// ** --> the strides of arrays can change when accessed by different routines for different purposes

// N_SLOPE_PPM          : size of Slope_PPM[]
// N_FC_VAR             : size of FC_Var[]
// N_FC_FLUX            : size of FC_Flux[]
// N_FL_FLUX/N_HF_FLUX  : for accessing FC_Flux[]
//                        --> may be different from N_FC_FLUX
//                            --> for example, in MHM_RP FC_Flux[] is also linked to Half_Flux[] used by
//                                Hydro_RiemannPredict_Flux() and Hydro_RiemannPredict()
//                            --> for the latter two routines, Half_Flux[] is accessed with N_HF_FLUX
//                                that is smaller than N_FC_FLUX
// N_HF_VAR             : for accessing PriVar_Half[], which is linked to PriVar[] with the size FLU_NXT^3
//                        --> also for accessing FC_B_Half[] in MHD

// NWAVE                : number of characteristic waves
// NCOMP_TOTAL_PLUS_MAG : total number of fluid variables plus magnetic field
// MAG_OFFSET           : array offset of magnetic field for arrays with the size NCOMP_TOTAL_PLUS_MAG
//
#define NCOMP_TOTAL_PLUS_MAG     ( NCOMP_TOTAL + NCOMP_MAG )

#ifdef MHD
#  define NWAVE                  ( NCOMP_FLUID + 2 )
#  define MAG_OFFSET             ( NCOMP_TOTAL )
#else
#  define NWAVE                  ( NCOMP_FLUID )
#  define MAG_OFFSET             ( NULL_INT )
#endif

#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  if   ( FLU_SCHEME == MHM )

#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 1 )
#     define N_FC_FLUX           ( N_FL_FLUX )

#  elif ( FLU_SCHEME == MHM_RP )

#    ifdef MHD
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 2 )
#     define N_HF_FLUX           ( FLU_NXT )
#    else
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( PS2 + 1 )
#     define N_HF_FLUX           ( FLU_NXT - 1 )
#    endif
#     define N_FC_FLUX           ( N_HF_FLUX )
#     define N_HF_VAR            ( FLU_NXT - 2 )

#  elif ( FLU_SCHEME == CTU )

#    ifdef MHD
#     define N_FC_VAR            ( PS2 + 4 )
#     define N_FL_FLUX           ( N_FC_VAR - 2 )
#     define N_HF_VAR            ( PS2 + 2 )
#    else
#     define N_FC_VAR            ( PS2 + 2 )
#     define N_FL_FLUX           ( N_FC_VAR )
#    endif
#     define N_HF_FLUX           ( N_FC_VAR )
#     define N_FC_FLUX           ( N_HF_FLUX )

#  endif // FLU_SCHEME

#  define N_SLOPE_PPM            ( N_FC_VAR + 2 )

#  ifdef MHD
#   define N_HF_ELE              ( N_FC_FLUX - 1 )
#   define N_FL_ELE              ( N_FL_FLUX - 1 )
#   define N_EC_ELE              ( N_FC_FLUX - 1 )
#  else
#   define N_EC_ELE              0
#  endif

#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )


// check non-physical negative values (e.g., negative density) for the fluid solver
#if ( defined GAMER_DEBUG  &&  MODEL == HYDRO )
#  define CHECK_UNPHYSICAL_IN_FLUID
#endif

#ifdef CHECK_UNPHYSICAL_IN_FLUID
#  include "stdio.h"
   bool Hydro_CheckUnphysical( const real Cons[], const real Prim[], const real* const Input, const real Passive[],
                               const char Info[], const char File[], const char Function[], const int Line, bool Show );
#endif


// perform spatial data reconstruction in characteristic variables (default: primitive variables)
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
//#  define CHAR_RECONSTRUCTION
#endif


// perform spatial data reconstruction in internal energy and use that to convert the face-centered
// primitive variables to conservative variables
// --> when it's disabled, the internal energy is converted from pressure using a given EoS
// --> pros: improve performance since the EoS conversion could be expensive
//     cons: (1) face-centered internal energy and pressure will NOT be fully self-consistent for a given EoS,
//               which may affect the accuracy of Riemann solver if it loads both conservative and primitive variables
//           (2) incompatible with CTU as it requires applying the characteristic tracing step to internal energy,
//               which has not been implemented
// --> unnecessary for EOS_GAMMA/EOS_ISOTHERMAL as they are fast
// --> disable it by default
#if ( EOS != EOS_GAMMA  &&  EOS != EOS_ISOTHERMAL  &&  FLU_SCHEME != CTU )
//#  define LR_EINT
#endif


// total number of target variables in the data reconstruction
#ifdef LR_EINT
#  define NCOMP_LR   ( NCOMP_TOTAL_PLUS_MAG + 1 )
#else
#  define NCOMP_LR   ( NCOMP_TOTAL_PLUS_MAG     )
#endif


// verify that the density and pressure in the intermediate states of Roe's Riemann solver are positive.
// --> if either is negative, we switch to other Riemann solvers (EXACT/HLLE/HLLC/HLLD)
#if (  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  &&  RSOLVER == ROE  )
#  ifdef MHD
//#     define CHECK_INTERMEDIATE    HLLD
#     define CHECK_INTERMEDIATE    HLLE
#  else
//#     define CHECK_INTERMEDIATE    HLLC
#     define CHECK_INTERMEDIATE    HLLE
#  endif
#endif


// use Eulerian with Y factor for Roe Solver in MHD
#if (  defined MHD  &&  ( RSOLVER == ROE || RSOLVER == HLLE )  )
#  define EULERY
#endif


// do not use the reference states for HLL solvers during the data reconstruction, as suggested in ATHENA
#if (  defined RSOLVER  &&  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  )

#  define HLL_NO_REF_STATE

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  ifdef HLL_NO_REF_STATE
#     define HLL_INCLUDE_ALL_WAVES
#  endif

#endif


// wave-speed estimates in the HLL-like Riemann solvers
#define HLL_WAVESPEED_ROE     1  // Roe average eigenvalues (Batten et al. 1997, SIAM J. Sci. Comput., 18, 1553)
#define HLL_WAVESPEED_PVRS    2  // Primitive Variable Riemann Solver (Toro 1999, Sec. 10.5.2)
#define HLL_WAVESPEED_DAVIS   3  // min/max of the left and right eigenvalues (Davis 1988, SIAM J. Sci. Stat, Comput., 9, 445)

// supported options:
// -> HLL_WAVESPEED_ROE (1) only supports the constant-gamma EoS (i.e., EOS_GAMMA)
//    HLL_WAVESPEED_PVRS (2) does not support MHD
// -> HLLC:
//       MHD on : none
//       MHD off: 1 for constant-gamma EoS and 2/3 for all EoS
//    HLLE:
//       MHD on : 1 for constant-gamma EoS and 3 for all EoS
//       MHD off: 1 for constant-gamma EoS and 2/3 for all EoS
//    HLLD:
//       MHD on : 3 for all EoS
//       MHD off: none

#  define HLLC_WAVESPEED   HLL_WAVESPEED_DAVIS
//#  define HLLC_WAVESPEED   HLL_WAVESPEED_PVRS
#ifdef MHD
#  define HLLE_WAVESPEED   HLL_WAVESPEED_DAVIS
#else
#  define HLLE_WAVESPEED   HLL_WAVESPEED_DAVIS
//#  define HLLE_WAVESPEED   HLL_WAVESPEED_PVRS
#endif
#  define HLLD_WAVESPEED   HLL_WAVESPEED_DAVIS


// 2. ELBDM macro
//=========================================================================================
#elif ( MODEL == ELBDM )


#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL



// ###############################################################################################
// number of threads in x/y directions for different solvers:
// Better to be a multiple of myWarpAllocationGranularity (WAG) * limitThreadsPerWarp(TPW)
// --> Please refer to CUDA_Occupancy_Calculator
//
// ComputeCompability  WAG   TPW   WAG*TPW
//       2.x             2    32        64
//       3.x             4    32       128
//       5.x             4    32       128
// ###############################################################################################

// 1. hydro solver
//=========================================================================================
#if ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )

#     define FLU_BLOCK_SIZE_X       FLU_NXT

#  ifdef FLOAT8
#     define FLU_BLOCK_SIZE_Y       4
#  else
#     define FLU_BLOCK_SIZE_Y       8
#  endif

#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == AMPERE )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  else
#     define FLU_BLOCK_SIZE_X       NULL_INT
#     ifdef GPU
#     error : UNKNOWN GPU_ARCH !!
#     endif
#  endif

#     define FLU_BLOCK_SIZE_Y       1

#elif ( FLU_SCHEME == CTU )

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512
#     endif
#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  elif ( GPU_ARCH == AMPERE )
#     ifdef FLOAT8
#     define FLU_BLOCK_SIZE_X       256
#     else
#     define FLU_BLOCK_SIZE_X       512      // not optimized yet
#     endif
#  else
#     define FLU_BLOCK_SIZE_X       NULL_INT
#     ifdef GPU
#     error : UNKNOWN GPU_ARCH !!
#     endif
#  endif

#     define FLU_BLOCK_SIZE_Y       1

#else
#  error : ERROR : unsupported hydro scheme in the makefile !!
#endif


// 2. ELBDM kinematic solver
//=========================================================================================
#elif ( MODEL == ELBDM )
#     define FLU_BLOCK_SIZE_X       PS2

#  if   ( GPU_ARCH == FERMI )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    4
#     else
#        define FLU_BLOCK_SIZE_Y    8
#     endif

#  elif ( GPU_ARCH == KEPLER )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == MAXWELL )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == PASCAL )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == VOLTA )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == TURING )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  elif ( GPU_ARCH == AMPERE )
#     ifdef FLOAT8
#        define FLU_BLOCK_SIZE_Y    16    // not optimized yet
#     else
#        define FLU_BLOCK_SIZE_Y    32    // not optimized yet
#     endif

#  else
#        define FLU_BLOCK_SIZE_Y    NULL_INT
#        ifdef GPU
#        error : UNKNOWN GPU_ARCH !!
#        endif
#  endif

#endif // MODEL


// 3. dt solver for fluid
//=========================================================================================
#     define DT_FLU_BLOCK_SIZE      512

// use shuffle reduction in the KEPLER and later GPUs
#  if ( GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA  ||  \
        GPU_ARCH == TURING  ||  GPU_ARCH == AMPERE )
#     define DT_FLU_USE_SHUFFLE
#  endif


// 4. source-term solver
//=========================================================================================
#     define SRC_BLOCK_SIZE         512



// warp size (which must be the same as the CUDA predefined constant "warpSize")
// --> please refer to https://en.wikipedia.org/wiki/CUDA#Version_features_and_specifications
//     for information on warp size
#ifdef __CUDACC__
#if ( GPU_ARCH == FERMI  ||  GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA  ||  \
      GPU_ARCH == TURING  ||  GPU_ARCH == AMPERE  )
// CUPOT.h will define WARP_SIZE as well
#  ifndef WARP_SIZE
#  define WARP_SIZE 32
#  endif
#elif defined GPU
#  error : UNKNOWN GPU_ARCH !!
#endif
#endif // #ifdef __CUDACC__



// #########################
// ## CPU/GPU integration ##
// #########################

// GPU device function specifier
#ifdef __CUDACC__
# define GPU_DEVICE          __forceinline__ __device__
# define GPU_DEVICE_NOINLINE    __noinline__ __device__
#else
# define GPU_DEVICE
# define GPU_DEVICE_NOINLINE
#endif

// unified CPU/GPU loop
#ifdef __CUDACC__
# define CGPU_LOOP( var, niter )    for (int (var)=threadIdx.x; (var)<(niter); (var)+=blockDim.x)
#else
# define CGPU_LOOP( var, niter )    for (int (var)=0;           (var)<(niter); (var)++          )
#endif



#endif // #ifndef __CUFLU_H__
