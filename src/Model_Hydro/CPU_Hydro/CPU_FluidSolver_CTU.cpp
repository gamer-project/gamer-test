#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  FLU_SCHEME == CTU && !defined SRHD )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#ifdef MHD
#include "CUFLU_Shared_ConstrainedTransport.cu"
#endif

#include "CUDA_ConstMemory.h"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh,
                               const real MinDens, const real MinPres, const real MinEint,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff,
                               const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2E_t EoS_DensPres2Eint,
                               const EoS_DP2C_t EoS_DensPres2CSqr,
                               const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const EoS_TEM2H_t EoS_Temp2HTilde,
                               const double EoS_AuxArray_Flt[],
                               const int    EoS_AuxArray_Int[],
                               const real *const EoS_Table[EOS_NTABLE_MAX] );
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time, const bool UsePot,
                        const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func, const double ExtAcc_AuxArray[],
                        const real MinDens, const real MinPres, const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ],
                        const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2C_t EoS_DensPres2CSqr,
                        const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp, EoS_TEM2H_t EoS_Temp2HTilde,
                        const EoS_TEM2C_t EoS_Temper2CSqr, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                        const real *const EoS_Table[EOS_NTABLE_MAX] );
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_FC_B[][ PS2P1*SQR(PS2) ], const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real dt, const real dh, const real MinDens, const real MinEint,
                           const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                           const double EoS_AuxArray_Flt[] );
#ifdef MHD
void MHD_ComputeElectric(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                          const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                          const real g_PriVar[][ CUBE(FLU_NXT) ],
                          const int NEle, const int NFlux, const int NPri, const int OffsetPri,
                          const real dt, const real dh,
                          const bool DumpIntEle, real g_IntEle[][NCOMP_ELE][ PS2P1*PS2 ],
                          const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[], const double Time,
                          const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                          const double ExtAcc_AuxArray[] );
void MHD_UpdateMagnetic( real *g_FC_Bx_Out, real *g_FC_By_Out, real *g_FC_Bz_Out,
                         const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                         const real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                         const real dt, const real dh, const int NOut, const int NEle, const int Offset_B_In );
void MHD_HalfStepPrimitive( const real g_Flu_In[][ CUBE(FLU_NXT) ],
                            const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                  real g_PriVar_Out[][ CUBE(FLU_NXT) ],
                            const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real dt, const real dh, const real MinDens );
#endif // #ifdef MHD

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
GPU_DEVICE
void Hydro_TGradientCorrection(       real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)  ],
                                const real g_FC_Flux  [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B_In  [][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                                const real g_PriVar   [][ CUBE(FLU_NXT) ],
                                const real dt, const real dh, const real MinDens, const real MinEint );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_FluidSolver_CTU
// Description :  CPU/GPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  1. Ref: (a) Stone et al., ApJS, 178, 137 (2008)
//                        (b) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                2. See include/CUFLU.h for the values and description of different symbolic constants
//                   such as N_FC_VAR, N_FC_FLUX, N_SLOPE_PPM, N_FL_FLUX, N_HF_VAR
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_In         : Array storing the input fluid variables
//                g_Flu_Array_Out        : Array to store the output fluid variables
//                g_Mag_Array_In         : Array storing the input B field (for MHD only)
//                g_Mag_Array_Out        : Array to store the output B field (for MHD only)
//                g_DE_Array_Out         : Array to store the dual-energy status
//                g_Flux_Array           : Array to store the output fluxes
//                g_Ele_Array            : Array to store the output electric field (for MHD only)
//                g_Corner_Array         : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                g_Pot_Array_USG        : Array storing the input potential for UNSPLIT_GRAVITY
//                g_PriVar               : Array to store the primitive variables
//                g_Slope_PPM            : Array to store the slope for the PPM reconstruction
//                g_FC_Var               : Array to store the half-step variables
//                g_FC_Flux              : Array to store the face-centered fluxes
//                g_FC_Mag_Half          : Array to store the half-step B field (for MHD only)
//                g_EC_Ele               : Array to store the edge-centered electric field (for MHD only)
//                NPatchGroup            : Number of patch groups to be evaluated
//                dt                     : Time interval to advance solution
//                dh                     : Cell size
//                StoreFlux              : true --> store the coarse-fine fluxes
//                StoreElectric          : true --> store the coarse-fine electric field
//                LR_Limiter             : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                         (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                        vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff           : Coefficient of the generalized MinMod limiter
//                Time                   : Current physical time                                 (for UNSPLIT_GRAVITY only)
//                UsePot                 : Add self-gravity and/or external potential            (for UNSPLIT_GRAVITY only)
//                ExtAcc                 : Add external acceleration                             (for UNSPLIT_GRAVITY only)
//                ExtAcc_Func            : Function pointer to the external acceleration routine (for UNSPLIT_GRAVITY only)
//                c_ExtAcc_AuxArray      : Auxiliary array for adding external acceleration      (for UNSPLIT_GRAVITY and CPU only)
//                                         --> When using GPU, this array is stored in the constant memory header
//                                             CUDA_ConstMemory.h and does not need to be passed as a function argument
//                MinDens/Pres/Eint      : Density, pressure, and internal energy floors
//                DualEnergySwitch       : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive            : true --> normalize passive scalars so that the sum of their mass density
//                                                  is equal to the gas mass density
//                NNorm                  : Number of passive scalars to be normalized
//                                         --> Should be set to the global variable "PassiveNorm_NVar"
//                c_NormIdx              : Target variable indices to be normalized
//                                         --> Should be set to the global variable "PassiveNorm_VarIdx"
//                                         --> When using GPU, this array is stored in the constant memory and does
//                                             not need to be passed as a function argument
//                                             --> Declared in CUDA_ConstMemory.h with the prefix "c_" to
//                                                 highlight that this is a constant variable on GPU
//                JeansMinPres           : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff     : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS_DensEint2Pres_Func : Function pointer to the EoS routine of computing the gas pressure
//                EoS_DensPres2Eint_Func :                    . . .                             gas internal energy
//                EoS_DensPres2CSqr_Func :                    . . .                             sound speed square
//                c_EoS_AuxArray_*       : Auxiliary arrays for the EoS routines (for CPU only)
//                c_EoS_Table            : EoS tables                            (for CPU only)
//                                         --> When using GPU, these CPU-only variables are stored in the constant memory
//                                             header CUDA_ConstMemory.h and do not need to be passed as function arguments
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
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
   const EoS_DP2C_t EoS_DensPres2CSqr_Func )
#else
void CPU_FluidSolver_CTU(
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
   const int NPatchGroup,
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm,
   const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func,
   const double c_EoS_AuxArray_Flt[],
   const int    c_EoS_AuxArray_Int[],
   const real* const c_EoS_Table[EOS_NTABLE_MAX] )
#endif // #ifdef __CUDACC__ ... else ...
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel          = true;
#  else
   const bool CorrHalfVel          = false;
#  endif
   const bool CorrHalfVel_No       = false;
   const bool StoreFlux_No         = false;
   const bool Con2Pri_Yes          = true;
#  ifdef MHD
   const bool StoreElectric_No     = false;
#  endif
#  if ( defined __CUDACC__  &&  !defined GRAVITY )
   const double *c_ExtAcc_AuxArray = NULL;
#  endif


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//    point to the arrays associated with different OpenMP threads (for CPU) or CUDA thread blocks (for GPU)
#     ifdef __CUDACC__
      const int array_idx = blockIdx.x;
#     else
#     ifdef OPENMP
      const int array_idx = omp_get_thread_num();
#     else
      const int array_idx = 0;
#     endif
#     endif // #ifdef __CUDACC__ ... else ...

      real (*const g_FC_Var_1PG   )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ] = g_FC_Var   [array_idx];
      real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [array_idx];
      real (*const g_PriVar_1PG   )                      [ CUBE(FLU_NXT)     ] = g_PriVar   [array_idx];
      real (*const g_Slope_PPM_1PG)[NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[array_idx];

#     ifdef MHD
      real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = g_FC_Mag_Half[array_idx];
      real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = g_EC_Ele     [array_idx];
#     else
      real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
      real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = NULL;
#     endif

#     ifdef MHD
      real (*const g_PriVar_Half_1PG)[ CUBE(FLU_NXT) ] = g_PriVar_1PG;
#     endif


//    loop over all patch groups
//    --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//        to work on different patch groups
#     ifdef __CUDACC__
      const int P = blockIdx.x;
#     else
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
#     endif
      {
//       1. evaluate the face-centered values at the half time-step
         Hydro_DataReconstruction( g_Flu_Array_In[P], g_Mag_Array_In[P], g_PriVar_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                   Con2Pri_Yes, LR_Limiter, MinMod_Coeff, dt, dh,
                                   MinDens, MinPres, MinEint, NormPassive, NNorm, c_NormIdx,
                                   JeansMinPres, JeansMinPres_Coeff,
                                   EoS_DensEint2Pres_Func, EoS_DensPres2Eint_Func, EoS_DensPres2CSqr_Func,
                                   NULL, NULL, NULL, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );


//       2. evaluate the face-centered half-step fluxes by solving the Riemann problem
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_HF_FLUX, 0, 0, CorrHalfVel_No,
                            NULL, NULL, NULL_REAL, NULL_REAL, NULL_REAL,
                            EXT_POT_NONE, EXT_ACC_NONE, NULL, NULL,
                            MinDens, MinPres, StoreFlux_No, NULL,
                            EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func,
                            NULL, NULL, NULL, NULL, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );


//       3. evaluate electric field and update B field at the half time-step
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_1PG, N_HF_ELE, N_HF_FLUX,
                              FLU_NXT, LR_GHOST_SIZE, dt, dh, StoreElectric_No, NULL,
                              CorrHalfVel_No, NULL, NULL, NULL_REAL,
                              EXT_POT_NONE, EXT_ACC_NONE, NULL, NULL );

         MHD_UpdateMagnetic( g_FC_Mag_Half_1PG[0], g_FC_Mag_Half_1PG[1], g_FC_Mag_Half_1PG[2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, (real)0.5*dt, dh, N_HF_VAR, N_HF_ELE, FLU_GHOST_SIZE-1 );
#        endif


//       4. correct the face-centered variables by the transverse flux gradients
         Hydro_TGradientCorrection( g_FC_Var_1PG, g_FC_Flux_1PG, g_Mag_Array_In[P], g_FC_Mag_Half_1PG, g_EC_Ele_1PG, g_PriVar_1PG,
                                    dt, dh, MinDens, MinEint );


//       5. evaluate the cell-centered primitive variables at the half time-step
//          --> for computing CT electric field later
#        ifdef MHD
         MHD_HalfStepPrimitive( g_Flu_Array_In[P], g_FC_Mag_Half_1PG, g_PriVar_Half_1PG, g_FC_Flux_1PG, dt, dh, MinDens );
#        endif


//       6. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
#        ifdef MHD
         const int NSkip_N = 1;
         const int NSkip_T = 1;
#        else
         const int NSkip_N = 0;
         const int NSkip_T = 1;
#        endif
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_FL_FLUX, NSkip_N, NSkip_T, CorrHalfVel,
                            g_Pot_Array_USG[P], g_Corner_Array[P], dt, dh, Time,
                            UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray,
                            MinDens, MinPres, StoreFlux, g_Flux_Array[P],
                            EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func,
                            NULL, NULL, NULL, NULL,
                            c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );


//       7. evaluate electric field and update B field at the full time-step
//          --> must update B field before Hydro_FullStepUpdate() since the latter requires
//              the updated magnetic energy when adopting the dual-energy formalism
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_Half_1PG, N_FL_ELE, N_FL_FLUX,
                              N_HF_VAR, 0, dt, dh, StoreElectric, g_Ele_Array[P],
                              CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P], Time,
                              UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray );

         MHD_UpdateMagnetic( g_Mag_Array_Out[P][0], g_Mag_Array_Out[P][1], g_Mag_Array_Out[P][2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, dt, dh, PS2, N_FL_ELE, FLU_GHOST_SIZE );
#        endif


//       8. full-step evolution of the fluid data
         Hydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], g_DE_Array_Out[P], g_Mag_Array_Out[P],
                               g_FC_Flux_1PG, dt, dh, MinDens, MinEint, DualEnergySwitch,
                               NormPassive, NNorm, c_NormIdx, c_EoS_AuxArray_Flt );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_TGradientCorrection
// Description :  Correct the face-centered variables by the transverse flux gradients
//
// Note        :  1. Ref: (a) Stone et al., ApJS, 178, 137 (2008)
//                        (b) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                2. Assuming "N_FC_VAR == N_HF_FLUX"
//
// Parameter   :  g_FC_Var     : Array to store the input and output face-centered conserved variables
//                               --> Accessed with the stride N_FC_VAR
//                g_FC_Flux    : Array storing the input face-centered fluxes
//                               --> Accessed with the stride N_HF_FLUX
//                g_FC_B_In    : Array storing the input initial   face-centered B field
//                g_FC_B_Half  : Array storing the input half-step face-centered B field
//                g_EC_Ele     : Array storing the input edge-centered electric field
//                g_PriVar     : Array storing the input cell-centered primitive variables
//                dt           : Time interval to advance solution
//                dh           : Cell size
//                MinDens/Eint : Density and internal energy floors
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_TGradientCorrection(       real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)  ],
                                const real g_FC_Flux  [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B_In  [][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                                const real g_PriVar   [][ CUBE(FLU_NXT) ],
                                const real dt, const real dh, const real MinDens, const real MinEint )
{

   const int  didx_flux[3]   = { 1, N_HF_FLUX, SQR(N_HF_FLUX) };
   const real dt_dh2         = (real)0.5*dt/dh;
#  ifdef MHD
   const int  didx_b_in  [3] = { 1, FLU_NXT,  SQR(FLU_NXT)  };
   const int  didx_b_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  didx_ele   [3] = { 1, N_HF_ELE, SQR(N_HF_ELE) };
   const real _dh            = (real)1.0/dh;
   const real dt_dh4         = (real)0.25*dt*_dh;
   const real dt_2           = (real)0.5*dt;

   real PriVar_1Cell[NCOMP_FLUID+NCOMP_MAG], B_Face[NCOMP_MAG][2];   // [2]=left/right faces
#  endif

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      real fc_var[2][NCOMP_TOTAL_PLUS_MAG];  // [2]=left/right faces

#     ifdef MHD
      const int nskip[3] = { 1, 1, 1 };
#     else
      int nskip[3];
      switch ( d )
      {
         case 0 : nskip[0] = 0;  nskip[1] = 1;  nskip[2] = 1;  break;
         case 1 : nskip[0] = 1;  nskip[1] = 0;  nskip[2] = 1;  break;
         case 2 : nskip[0] = 1;  nskip[1] = 1;  nskip[2] = 0;  break;
      }
#     endif

      const int size_i  = ( N_FC_VAR - 2*nskip[0] );
      const int size_j  = ( N_FC_VAR - 2*nskip[1] );
      const int size_k  = ( N_FC_VAR - 2*nskip[2] );
      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx0, size_i*size_j*size_k )
      {
//       i/j/k0 start from zero
         const int i0         = idx0 % size_i;
         const int j0         = idx0 % size_ij / size_i;
         const int k0         = idx0 / size_ij;

         const int i_fc_var   = i0 + nskip[0];
         const int j_fc_var   = j0 + nskip[1];
         const int k_fc_var   = k0 + nskip[2];
         const int idx_fc_var = IDX321( i_fc_var, j_fc_var, k_fc_var, N_FC_VAR, N_FC_VAR );

         const int idx_fluxR  = idx_fc_var;  // assuming N_FC_VAR == N_HF_FLUX
         const int idx_fluxL1 = idx_fluxR - didx_flux[TDir1];
         const int idx_fluxL2 = idx_fluxR - didx_flux[TDir2];

//       0. load g_FC_Var[] to the local variable fc[] to reduce the GPU global memory access
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            fc_var[0][v] = g_FC_Var[faceL][v][idx_fc_var];
            fc_var[1][v] = g_FC_Var[faceR][v][idx_fc_var];
         }


//       1. calculate the transverse fluid flux gradients and update the corresponding face-centered fluid variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            real Correct, TGrad1, TGrad2;

            TGrad1  = g_FC_Flux[TDir1][v][idx_fluxR] - g_FC_Flux[TDir1][v][idx_fluxL1];
            TGrad2  = g_FC_Flux[TDir2][v][idx_fluxR] - g_FC_Flux[TDir2][v][idx_fluxL2];
            Correct = -dt_dh2*( TGrad1 + TGrad2 );

            fc_var[0][v] += Correct;
            fc_var[1][v] += Correct;
         }


#        ifdef MHD
//       2. correct the transverse B field
         const int idx_ele = IDX321( i0, j0, k0, N_HF_ELE, N_HF_ELE );

         for (int v=1; v<NCOMP_MAG; v++)
         {
            real Correct, TGrad1, TGrad2, Sign;

            const int TD1 = (d+v)%3;            // transverse direction 1
            const int TD2 = (d+2*v)%3;          // transverse direction 2
            const int TB  = TD1 + MAG_OFFSET;   // target transverse B field

            TGrad1  = g_EC_Ele[d][ idx_ele + didx_ele[TD2]                 ] - g_EC_Ele[d][ idx_ele                 ];
            TGrad2  = g_EC_Ele[d][ idx_ele + didx_ele[TD2] + didx_ele[TD1] ] - g_EC_Ele[d][ idx_ele + didx_ele[TD1] ];
            Sign    = (real)2.0*v - (real)3.0;  // v=1/2 --> sign=-1/+1
            Correct = Sign*dt_dh4*( TGrad2 + TGrad1 );

            fc_var[0][TB] += Correct;
            fc_var[1][TB] += Correct;
         } // for (int v=1; v<NCOMP_MAG; v++)


//       3. add the divergence(B) source terms
         Hydro_Rotate3D( fc_var[0], d, true, MAG_OFFSET );
         Hydro_Rotate3D( fc_var[1], d, true, MAG_OFFSET );

//       3-1. get the initial cell-centered primitive variables
         const int i_pri   = i_fc_var + LR_GHOST_SIZE;
         const int j_pri   = j_fc_var + LR_GHOST_SIZE;
         const int k_pri   = k_fc_var + LR_GHOST_SIZE;
         const int idx_pri = IDX321( i_pri, j_pri, k_pri, FLU_NXT, FLU_NXT );

//       skip passive scalars
         for (int v=0; v<NCOMP_FLUID; v++)   PriVar_1Cell[ v               ] = g_PriVar[ v              ][idx_pri];
         for (int v=0; v<NCOMP_MAG;   v++)   PriVar_1Cell[ v + NCOMP_FLUID ] = g_PriVar[ v + MAG_OFFSET ][idx_pri];

         Hydro_Rotate3D( PriVar_1Cell, d, true, NCOMP_FLUID );

//       3-2. get the initial face-centered B field
         const int idx_b_in[3] = { IDX321( i_pri, j_pri, k_pri, FLU_NXT_P1, FLU_NXT    ),
                                   IDX321( i_pri, j_pri, k_pri, FLU_NXT,    FLU_NXT_P1 ),
                                   IDX321( i_pri, j_pri, k_pri, FLU_NXT,    FLU_NXT    ) };

         B_Face[0][0] = g_FC_B_In[  d  ][ idx_b_in[  d  ]                    ];
         B_Face[0][1] = g_FC_B_In[  d  ][ idx_b_in[  d  ] + didx_b_in[  d  ] ];

         B_Face[1][0] = g_FC_B_In[TDir1][ idx_b_in[TDir1]                    ];
         B_Face[1][1] = g_FC_B_In[TDir1][ idx_b_in[TDir1] + didx_b_in[TDir1] ];

         B_Face[2][0] = g_FC_B_In[TDir2][ idx_b_in[TDir2]                    ];
         B_Face[2][1] = g_FC_B_In[TDir2][ idx_b_in[TDir2] + didx_b_in[TDir2] ];


//       3-3. add the divergence(B) source term
         const real Vy = PriVar_1Cell[ 2 ];
         const real Vz = PriVar_1Cell[ 3 ];
         const real Bx = PriVar_1Cell[ 0 + NCOMP_FLUID ];
         const real By = PriVar_1Cell[ 1 + NCOMP_FLUID ];
         const real Bz = PriVar_1Cell[ 2 + NCOMP_FLUID ];

         real dB[NCOMP_MAG], SrcFlu[NCOMP_FLUID-1], SrcMag[2], Vy_MinModBxz, Vz_MinModBxy;

         for (int v=0; v<NCOMP_MAG; v++)  dB[v] = ( B_Face[v][1] - B_Face[v][0] )*_dh;

#        define MINMOD( a , b )  (  ( (a)*(b)>(real)0.0 ) ? ( SIGN(a)*FMIN(FABS(a),FABS(b)) ) : (real)0.0  )
         Vy_MinModBxz = Vy*MINMOD( -dB[2], dB[0] );
         Vz_MinModBxy = Vz*MINMOD( -dB[1], dB[0] );
#        undef MINMOD

         SrcFlu[0] = dt_2*Bx*dB[0];
         SrcFlu[1] = dt_2*By*dB[0];
         SrcFlu[2] = dt_2*Bz*dB[0];
         SrcFlu[3] = dt_2*( By*Vy_MinModBxz + Bz*Vz_MinModBxy );

         SrcMag[0] = dt_2*Vy_MinModBxz;
         SrcMag[1] = dt_2*Vz_MinModBxy;

         for (int f=0; f<2; f++)
         {
            fc_var[f][ 1 ] += SrcFlu[0];
            fc_var[f][ 2 ] += SrcFlu[1];
            fc_var[f][ 3 ] += SrcFlu[2];
            fc_var[f][ 4 ] += SrcFlu[3];

            fc_var[f][ 1 + MAG_OFFSET ] += SrcMag[0];
            fc_var[f][ 2 + MAG_OFFSET ] += SrcMag[1];
         }


//       4. set the longitudinal B field to the half-step values updated by MHD_UpdateMagnetic()
         int idx_b_half;

         switch ( d )
         {
            case 0 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR+1, N_HF_VAR   );  break;
            case 1 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR,   N_HF_VAR+1 );  break;
            case 2 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR,   N_HF_VAR   );  break;
         }

         fc_var[0][MAG_OFFSET] = g_FC_B_Half[d][ idx_b_half                  ];
         fc_var[1][MAG_OFFSET] = g_FC_B_Half[d][ idx_b_half + didx_b_half[d] ];

         Hydro_Rotate3D( fc_var[0], d, false, MAG_OFFSET );
         Hydro_Rotate3D( fc_var[1], d, false, MAG_OFFSET );
#        endif // #ifdef MHD


//       5. apply density and internal energy floors
         for (int f=0; f<2; f++)
         {
#           ifdef MHD
            const real Bx   = fc_var[f][ MAG_OFFSET + 0 ];
            const real By   = fc_var[f][ MAG_OFFSET + 1 ];
            const real Bz   = fc_var[f][ MAG_OFFSET + 2 ];
            const real Emag= (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#           else
            const real Emag = NULL_REAL;
#           endif
            fc_var[f][0] = FMAX( fc_var[f][0], MinDens );
            fc_var[f][4] = Hydro_CheckMinEintInEngy( fc_var[f][0], fc_var[f][1], fc_var[f][2], fc_var[f][3], fc_var[f][4],
                                                     MinEint, Emag );
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
            fc_var[f][v] = FMAX( fc_var[f][v], TINY_NUMBER );
#           endif
         }

//       store the results to g_FC_Var[]
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            g_FC_Var[faceL][v][idx_fc_var] = fc_var[0][v];
            g_FC_Var[faceR][v][idx_fc_var] = fc_var[1][v];
         }

      } // CGPU_LOOP( idx0, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_TGradientCorrection



#endif // #if ( MODEL == HYDRO  &&  FLU_SCHEME == CTU )
