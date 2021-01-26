#include "CUAPI.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef GPU



// fluid solver prototypes in different models
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
__global__ void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                                         const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ],
                                         const real dh, const real Safety, const real MinPres,
                                         const EoS_DE2P_t EoS_DensEint2Pres_Func, const EoS_DP2C_t EoS_DensPres2CSqr_Func,
                                         const EoS_TEM2C_t EoS_Temper2CSqr_Func, const EoS_GUESS_t EoS_GuessHTilde_Func,
                                         const EoS_H2TEM_t EoS_HTilde2Temp_Func );
#ifdef GRAVITY
__global__
void CUPOT_dtSolver_HydroGravity( real g_dt_Array[], const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                  const double g_Corner_Array[][3],
                                  const real dh, const real Safety, const bool P5_Gradient,
                                  const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                                  const double ExtAcc_Time );
#endif

#elif ( MODEL == ELBDM )
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const real Taylor3_Coeff, const bool XYZ, const real MinDens );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL


#ifdef GRAVITY

// Poisson solver prototypes
#if   ( POT_SCHEME == SOR )
#ifdef USE_PSOLVER_10TO14
__global__ void CUPOT_PoissonSolver_SOR_10to14cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6,
                                                    const real Const, const IntScheme_t IntScheme );
#else
__global__ void CUPOT_PoissonSolver_SOR_16to18cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6,
                                                    const real Const, const IntScheme_t IntScheme );
#endif // #ifdef USE_PSOLVER_10TO14 ... else ...
#elif ( POT_SCHEME == MG )
__global__ void CUPOT_PoissonSolver_MG( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                        const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                              real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                        const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                                        const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                                        const IntScheme_t IntScheme );
#endif // POT_SCHEME


// Gravity solver prototypes in different models
#if   ( MODEL == HYDRO )
__global__
void CUPOT_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_Emag_Array   [][ CUBE(PS1) ],
   const real dt, const real dh, const bool P5_Gradient,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double TimeNew, const double TimeOld, const real MinEint,
   const EoS_GUESS_t EoS_GuessHTilde_Func,
   const EoS_TEM2H_t EoS_Temp2HTilde_Func,
   const EoS_H2TEM_t EoS_HTilde2Temp_Func );

#elif ( MODEL == ELBDM )
__global__
void CUPOT_ELBDMGravitySolver(       real g_Flu_Array[][GRA_NIN][ PS1*PS1*PS1 ],
                               const real g_Pot_Array[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                               const double g_Corner_Array[][3],
                               const real EtaDt, const real dh, const real Lambda );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

#endif // GRAVITY


// source-term solver prototype
__global__
void CUSRC_SrcSolver_IterateAllCells(
   const real g_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
         real g_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
   const real g_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
   const double g_Corner_Array[][3],
   const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
   const double TimeNew, const double TimeOld,
   const real MinDens, const real MinPres, const real MinEint,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Set_Default_GPU_Parameter
// Description :  Set several GPU parameters to the default values if they are not set in the input file
//
// Parameter   :  GPU_NStream     : Number of streams for the asynchronous memory copy in GPU
//                Flu_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the fluid solver
//                Pot_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the Poisson solver
//                Che_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the Grackle solver
//                Src_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the source-term solver
//-------------------------------------------------------------------------------------------------------
void CUAPI_Set_Default_GPU_Parameter( int &GPU_NStream, int &Flu_GPU_NPGroup, int &Pot_GPU_NPGroup, int &Che_GPU_NPGroup,
                                      int &Src_GPU_NPGroup )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// get the device ID
   int GetDeviceID = 999;
   CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );


// load the device properties
   cudaDeviceProp DeviceProp;
   CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );


// set the default GPU parameters
// (1) GPU_NSTREAM
   if ( GPU_NStream <= 0 )
   {
      if ( DeviceProp.deviceOverlap )
      {
#        if   ( MODEL == HYDRO )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 8;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == MAXWELL )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == PASCAL )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == VOLTA )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == TURING )
            GPU_NStream = 16;
#           else
#           error : UNKNOWN GPU_ARCH !!
#           endif

#        elif ( MODEL == ELBDM )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 8;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == MAXWELL )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == PASCAL )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == VOLTA )
            GPU_NStream = 16;
#           elif ( GPU_ARCH == TURING )
            GPU_NStream = 16;
#           else
#           error : ERROR : UNKNOWN GPU_ARCH !!
#           endif
#        else
#           error : ERROR : UNKNOWN MODEL !!
#        endif // MODEL
      } // if ( DeviceProp.deviceOverlap )

      else
         GPU_NStream = 1;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d"
                              " --> might be further fine-tuned\n", "GPU_NSTREAM", GPU_NSTREAM );
   } // if ( GPU_NStream <= 0 )


// (2) XXX_GPU_NPGROUP
// (2-1) FLU_GPU_NPGROUP
   if ( Flu_GPU_NPGroup <= 0 )
   {
#     if   ( MODEL == HYDRO )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == MAXWELL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == PASCAL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == VOLTA )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == TURING )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif

#     elif ( MODEL == ELBDM )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == MAXWELL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == PASCAL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == VOLTA )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == TURING )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif
#     else
#        error : ERROR : UNKNOWN MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d"
                              " --> might be further fine-tuned\n", "FLU_GPU_NPGROUP", Flu_GPU_NPGroup );
   } // if ( Flu_GPU_NPGroup <= 0 )

// (2-2) POT_GPU_NPGROUP
#  ifdef GRAVITY
   if ( Pot_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d"
                              " --> might be further fine-tuned\n", "POT_GPU_NPGROUP", Pot_GPU_NPGroup );
   } // if ( Pot_GPU_NPGroup <= 0 )
#  endif

// (2-3) CHE_GPU_NPGROUP
#  ifdef SUPPORT_GRACKLE
   if ( Che_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d"
                              " --> might be further fine-tuned\n", "CHE_GPU_NPGROUP", Che_GPU_NPGroup );
   } // if ( Che_GPU_NPGroup <= 0 )
#  endif

// (2-4) SRC_GPU_NPGROUP
   if ( Src_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d"
                              " --> might be further fine-tuned\n", "SRC_GPU_NPGROUP", Src_GPU_NPGroup );
   } // if ( Src_GPU_NPGroup <= 0 )


// (3) cache preference
// (3-1) fluid solver
#  if   ( MODEL == HYDRO )
#  if   ( FLU_SCHEME == RTVD )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_RTVD,             cudaFuncCachePreferShared )  );
#  elif ( FLU_SCHEME == MHM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,              cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == MHM_RP )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,              cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == CTU )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_CTU,              cudaFuncCachePreferL1     )  );
#  endif
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_dtSolver_HydroCFL,            cudaFuncCachePreferShared )  );
#  ifdef GRAVITY
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_dtSolver_HydroGravity,        cudaFuncCachePreferShared )  );
#  endif

#  elif ( MODEL == ELBDM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver,                  cudaFuncCachePreferShared )  );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY

// (3-2) Poisson solver
#  if   ( POT_SCHEME == SOR )
#  ifdef USE_PSOLVER_10TO14
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_PoissonSolver_SOR_10to14cube, cudaFuncCachePreferShared )  );
#  else
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_PoissonSolver_SOR_16to18cube, cudaFuncCachePreferShared )  );
#  endif
#  elif ( POT_SCHEME == MG )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_PoissonSolver_MG,             cudaFuncCachePreferShared )  );
#  endif // POT_SCHEME


// (3-3) gravity solver
#  if   ( MODEL == HYDRO )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_HydroGravitySolver,           cudaFuncCachePreferShared )  );

#  elif ( MODEL == ELBDM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_ELBDMGravitySolver,           cudaFuncCachePreferL1     )  );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  endif // GRAVITY


// (3-4) source-term solver
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUSRC_SrcSolver_IterateAllCells,   cudaFuncCachePreferL1      )  );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_Set_Default_GPU_Parameter



#endif // #ifdef GPU
