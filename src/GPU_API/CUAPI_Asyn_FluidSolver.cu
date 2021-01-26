#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
__global__ void CUFLU_FluidSolver_RTVD(
   real g_Fluid_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real g_Fluid_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real g_Flux     [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner[][3],
   const real g_Pot_USG[][ CUBE(USG_NXT_F) ],
   const real dt, const real _dh, const bool StoreFlux,
   const bool XYZ, const real MinDens, const real MinPres, const real MinEint,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
__global__
void CUFLU_FluidSolver_MHM(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func,
   const EoS_GUESS_t EoS_GuessHTilde_Func,
   const EoS_H2TEM_t EoS_HTilde2Temp_Func,
   const EoS_TEM2H_t EoS_Temp2HTilde_Func,
   const EoS_TEM2C_t EoS_Temper2CSqr_Func );
#elif ( FLU_SCHEME == CTU )
__global__
void CUFLU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func );
#endif // FLU_SCHEME

#elif ( MODEL == ELBDM )
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const real Taylor3_Coeff, const bool XYZ, const real MinDens );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

#ifndef GRAVITY
static ExtAcc_t GPUExtAcc_Ptr = NULL;
#endif


// device pointers
extern real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ];
extern double (*d_Corner_Array_F)[3];
#if ( MODEL == HYDRO )
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ CUBE(PS2) ];
#else
static char (*d_DE_Array_F_Out)[ CUBE(PS2) ] = NULL;
#endif
#ifdef MHD
extern real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)         ];
extern real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ];
#else
static real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
static real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)          ] = NULL;
static real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ]            = NULL;
#endif
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*d_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*d_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*d_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#else
static real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
static real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ] = NULL;
#endif // MHD
#endif // FLU_SCHEME
#endif // #if ( MODEL == HYDRO )

#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ];
#else
static real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ] = NULL;
#endif

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_FluidSolver
// Description :  1. MODEL == HYDRO : use GPU to solve the Euler equations by different schemes
//                                    --> invoke the kernel "CUFLU_FluidSolver_XXX"
//                2. MODEL == ELBDM : use GPU to solve the kinematic operator in the Schrodinger's equations
//                                    --> invoke the kernel "CUFLU_ELBDMSolver"
//
//                ***********************************************************
//                **                Asynchronous Function                  **
//                **                                                       **
//                **  will return before the execution in GPU is complete  **
//                ***********************************************************
//
// Note        :  1. Use streams for the asychronous memory copy between device and host
//                2. Prefix "d" : for pointers pointing to the "Device" memory space
//                   Prefix "h" : for pointers pointing to the "Host"   memory space
//                3. Use the input pamameter "XYZ" to control the order of update for dimensional-splitting
//                   method (currently only RTVD)
//                4. Currently five hydro schemes are supported :
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   3. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   4. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
// Parameter   :  h_Flu_Array_In      : Host array to store the input fluid variables
//                h_Flu_Array_Out     : Host array to store the output fluid variables
//                h_Mag_Array_In      : Host array storing the input B field (for MHD only)
//                h_Mag_Array_Out     : Host array to store the output B field (for MHD only)
//                h_DE_Array_Out      : Host array to store the dual-energy status
//                h_Flux_Array        : Host array to store the output fluxes
//                h_Ele_Array         : Host array to store the output electric field (for MHD only)
//                h_Corner_Array      : Host array storing the physical corner coordinates of each patch group
//                h_Pot_Array_USG     : Host array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup         : Number of patch groups evaluated simultaneously by GPU
//                dt                  : Time interval to advance solution
//                dh                  : Cell size
//                StoreFlux           : true --> store the coarse-fine fluxes
//                StoreElectric       : true --> store the coarse-fine electric field
//                XYZ                 : true  : x->y->z ( forward sweep)
//                                      false : z->y->x (backward sweep)
//                                      ~ useless in directionally unsplit schemes
//                LR_Limiter          : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                      (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                     vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff        : Coefficient of the generalized MinMod limiter
//                ELBDM_Eta           : Particle mass / Planck constant
//                ELBDM_Taylor3_Coeff : Coefficient in front of the third term in the Taylor expansion for ELBDM
//                ELBDM_Taylor3_Auto  : true --> Determine ELBDM_Taylor3_Coeff automatically by invoking the
//                                               function "ELBDM_SetTaylor3Coeff"
//                Time                : Current physical time                      (for UNSPLIT_GRAVITY only)
//                UsePot              : Add self-gravity and/or external potential (for UNSPLIT_GRAVITY only)
//                ExtAcc              : Add external acceleration                  (for UNSPLIT_GRAVITY only)
//                MinDens/Pres/Eint   : Density, pressure, and internal energy floors
//                DualEnergySwitch    : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive         : true --> normalize passive scalars so that the sum of their mass density
//                                               is equal to the gas mass density
//                NNorm               : Number of passive scalars to be normalized
//                                      --> Should be set to the global variable "PassiveNorm_NVar"
//                JeansMinPres        : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff  : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                GPU_NStream         : Number of CUDA streams for the asynchronous memory copy
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In[][FLU_NIN ][ CUBE(FLU_NXT) ],
                             real h_Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                             real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                             real h_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                             char h_DE_Array_Out[][ CUBE(PS2) ],
                             real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                             real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                             const double h_Corner_Array[][3],
                             real h_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
                             const int NPatchGroup, const real dt, const real dh,
                             const bool StoreFlux, const bool StoreElectric,
                             const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real ELBDM_Eta, real ELBDM_Taylor3_Coeff, const bool ELBDM_Taylor3_Auto,
                             const double Time, const bool UsePot, const OptExtAcc_t ExtAcc,
                             const real MinDens, const real MinPres, const real MinEint,
                             const real DualEnergySwitch, const bool NormPassive, const int NNorm,
                             const bool JeansMinPres, const real JeansMinPres_Coeff,
                             const int GPU_NStream )
{

// check
#  ifdef GAMER_DEBUG
#  if   ( MODEL == HYDRO )
   if ( LR_Limiter != VANLEER  &&  LR_Limiter != GMINMOD  &&  LR_Limiter != ALBADA  &&  LR_Limiter != EXTPRE  &&
        LR_Limiter != VL_GMINMOD  &&  LR_Limiter != LR_LIMITER_NONE )
      Aux_Error( ERROR_INFO, "unsupported limiter (%d) !!\n", LR_Limiter );

#  ifdef UNSPLIT_GRAVITY
   if ( UsePot )
   {
      if ( h_Pot_Array_USG   == NULL )   Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );
      if ( d_Pot_Array_USG_F == NULL )   Aux_Error( ERROR_INFO, "d_Pot_Array_USG_F == NULL !!\n" );
   }

   if ( ExtAcc )
   {
      if ( h_Corner_Array   == NULL )    Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
      if ( d_Corner_Array_F == NULL )    Aux_Error( ERROR_INFO, "d_Corner_Array_F == NULL !!\n" );
   }
#  endif

#  elif ( MODEL == ELBDM )

#  else
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif

   if ( StoreFlux )
   {
      if ( d_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "d_Flux_Array == NULL !!\n" );
      if ( h_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "h_Flux_Array == NULL !!\n" );
   }

#  ifdef MHD
   if ( h_Mag_Array_In    == NULL ) Aux_Error( ERROR_INFO, "h_Mag_Array_In == NULL !!\n" );
   if ( d_Mag_Array_F_In  == NULL ) Aux_Error( ERROR_INFO, "d_Mag_Array_F_In == NULL !!\n" );

   if ( h_Mag_Array_Out   == NULL ) Aux_Error( ERROR_INFO, "h_Mag_Array_Out == NULL !!\n" );
   if ( d_Mag_Array_F_Out == NULL ) Aux_Error( ERROR_INFO, "d_Mag_Array_F_Out == NULL !!\n" );

   if ( d_FC_Mag_Half     == NULL ) Aux_Error( ERROR_INFO, "d_FC_Mag_Half == NULL !!\n" );
   if ( d_EC_Ele          == NULL ) Aux_Error( ERROR_INFO, "d_EC_Ele == NULL !!\n" );

   if ( StoreElectric )
   {
      if ( d_Ele_Array == NULL )   Aux_Error( ERROR_INFO, "d_Ele_Array == NULL !!\n" );
      if ( h_Ele_Array == NULL )   Aux_Error( ERROR_INFO, "h_Ele_Array == NULL !!\n" );
   }
#  endif
#  endif // #ifdef GAMER_DEBUG


   const dim3 BlockDim_FluidSolver ( FLU_BLOCK_SIZE_X, FLU_BLOCK_SIZE_Y, 1 ); // for the fluidsolvers

// model-dependent operations
#  if   ( MODEL == HYDRO )

#  elif ( MODEL == ELBDM )
// evaluate the optimized Taylor expansion coefficient
   if ( ELBDM_Taylor3_Auto )  ELBDM_Taylor3_Coeff = ELBDM_SetTaylor3Coeff( dt, dh, ELBDM_Eta );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   int *NPatch_per_Stream  = new int [GPU_NStream];
   int *UsedPatch          = new int [GPU_NStream];
   int *Flu_MemSize_In     = new int [GPU_NStream];
   int *Flu_MemSize_Out    = new int [GPU_NStream];
   int *Flux_MemSize       = new int [GPU_NStream];
#  ifdef MHD
   int *Mag_MemSize_In     = new int [GPU_NStream];
   int *Mag_MemSize_Out    = new int [GPU_NStream];
   int *Ele_MemSize        = new int [GPU_NStream];
#  endif
#  ifdef UNSPLIT_GRAVITY
   int *USG_MemSize        = new int [GPU_NStream];
   int *Corner_MemSize     = new int [GPU_NStream];
#  endif
#  ifdef DUAL_ENERGY
   int *DE_MemSize_Out     = new int [GPU_NStream];
#  endif


// set the number of patches of each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatchGroup;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatchGroup / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatchGroup - UsedPatch[GPU_NStream-1];
   }


// set the size of data to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Flu_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN *CUBE(FLU_NXT);
      Flu_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NOUT*CUBE(PS2);
      Flux_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*NFLUX_TOTAL*9*PS2*PS2;
#     ifdef MHD
      Mag_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_MAG*FLU_NXT_P1*SQR(FLU_NXT);
      Mag_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_MAG*PS2P1*SQR(PS2);
      Ele_MemSize    [s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_ELE*9*PS2P1*PS2;
#     endif
#     ifdef UNSPLIT_GRAVITY
      USG_MemSize    [s] = sizeof(real  )*NPatch_per_Stream[s]*CUBE(USG_NXT_F);
      Corner_MemSize [s] = sizeof(double)*NPatch_per_Stream[s]*3;
#     endif
#     ifdef DUAL_ENERGY
      DE_MemSize_Out [s] = sizeof(char  )*NPatch_per_Stream[s]*CUBE(PS2);
#     endif
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_F_In  + UsedPatch[s], h_Flu_Array_In  + UsedPatch[s],
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );

#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Mag_Array_F_In  + UsedPatch[s], h_Mag_Array_In  + UsedPatch[s],
                         Mag_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif

#     ifdef UNSPLIT_GRAVITY
      if ( UsePot )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_USG_F + UsedPatch[s], h_Pot_Array_USG + UsedPatch[s],
                         USG_MemSize   [s], cudaMemcpyHostToDevice, Stream[s] )  );

      if ( ExtAcc )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_F  + UsedPatch[s], h_Corner_Array  + UsedPatch[s],
                         Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if   ( MODEL == HYDRO )

#        if   ( FLU_SCHEME == RTVD )

         CUFLU_FluidSolver_RTVD <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              dt, 1.0/dh, StoreFlux, XYZ, MinDens, MinPres, MinEint,
              EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

         CUFLU_FluidSolver_MHM <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Mag_Array_F_In  + UsedPatch[s],
              d_Mag_Array_F_Out + UsedPatch[s],
              d_DE_Array_F_Out  + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Ele_Array       + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              d_PriVar          + UsedPatch[s],
              d_Slope_PPM       + UsedPatch[s],
              d_FC_Var          + UsedPatch[s],
              d_FC_Flux         + UsedPatch[s],
              d_FC_Mag_Half     + UsedPatch[s],
              d_EC_Ele          + UsedPatch[s],
              dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff,
              Time, UsePot, ExtAcc, GPUExtAcc_Ptr, MinDens, MinPres, MinEint,
              DualEnergySwitch, NormPassive, NNorm, JeansMinPres, JeansMinPres_Coeff,
              EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr,
              EoS_GuessHTilde_GPUPtr, EoS_HTilde2Temp_GPUPtr, EoS_Temp2HTilde_GPUPtr, EoS_Temper2CSqr_GPUPtr );

#        elif ( FLU_SCHEME == CTU )

         CUFLU_FluidSolver_CTU <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Mag_Array_F_In  + UsedPatch[s],
              d_Mag_Array_F_Out + UsedPatch[s],
              d_DE_Array_F_Out  + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Ele_Array       + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              d_PriVar          + UsedPatch[s],
              d_Slope_PPM       + UsedPatch[s],
              d_FC_Var          + UsedPatch[s],
              d_FC_Flux         + UsedPatch[s],
              d_FC_Mag_Half     + UsedPatch[s],
              d_EC_Ele          + UsedPatch[s],
              dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff,
              Time, UsePot, ExtAcc, GPUExtAcc_Ptr, MinDens, MinPres, MinEint,
              DualEnergySwitch, NormPassive, NNorm, JeansMinPres, JeansMinPres_Coeff,
              EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );

#        else

#        error : unsupported GPU hydro scheme

#        endif // FLU_SCHEME

#     elif ( MODEL == ELBDM )

         CUFLU_ELBDMSolver <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              dt, 1.0/dh, ELBDM_Eta, StoreFlux, ELBDM_Taylor3_Coeff, XYZ, MinDens );

#     else

#        error : unsupported MODEL !!

#     endif // MODEL


      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array_Out + UsedPatch[s], d_Flu_Array_F_Out + UsedPatch[s],
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( StoreFlux )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flux_Array    + UsedPatch[s], d_Flux_Array      + UsedPatch[s],
                         Flux_MemSize[s],    cudaMemcpyDeviceToHost, Stream[s] )  );

#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Mag_Array_Out + UsedPatch[s], d_Mag_Array_F_Out + UsedPatch[s],
                         Mag_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( StoreElectric )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Ele_Array     + UsedPatch[s], d_Ele_Array       + UsedPatch[s],
                         Ele_MemSize[s],    cudaMemcpyDeviceToHost, Stream[s] )  );
#     endif

#     ifdef DUAL_ENERGY
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_DE_Array_Out  + UsedPatch[s], d_DE_Array_F_Out  + UsedPatch[s],
                         DE_MemSize_Out[s],  cudaMemcpyDeviceToHost, Stream[s] )  );
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] Flu_MemSize_In;
   delete [] Flu_MemSize_Out;
   delete [] Flux_MemSize;
#  ifdef MHD
   delete [] Mag_MemSize_In;
   delete [] Mag_MemSize_Out;
   delete [] Ele_MemSize;
#  endif
#  ifdef UNSPLIT_GRAVITY
   delete [] USG_MemSize;
   delete [] Corner_MemSize;
#  endif
#  ifdef DUAL_ENERGY
   delete [] DE_MemSize_Out;
#  endif

} // FUNCTION : CUAPI_Asyn_FluidSolver



#endif // #ifdef GPU
