#include "CUAPI.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef GPU



#if   ( MODEL == HYDRO )
__global__
void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
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

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

// device pointers
extern real *d_dt_Array_T;
extern real (*d_Flu_Array_T)[FLU_NIN_T][ CUBE(PS1) ];
#ifdef GRAVITY
extern real (*d_Pot_Array_T)[ CUBE(GRA_NXT) ];
extern double (*d_Corner_Array_PGT)[3];
#endif
#ifdef MHD
extern real (*d_Mag_Array_T)[NCOMP_MAG][ PS1P1*SQR(PS1) ];
#else
static real (*d_Mag_Array_T)[NCOMP_MAG][ PS1P1*SQR(PS1) ] = NULL;
#endif

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_dtSolver
// Description :  Invoke various dt solvers
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
//
// Parameter   :  TSolver        : Target dt solver
//                                 --> DT_FLU_SOLVER : dt solver for fluid
//                                     DT_GRA_SOLVER : dt solver for gravity
//                h_dt_Array     : Host array to store the minimum dt in each target patch
//                h_Flu_Array    : Host array storing the prepared fluid     data of each target patch
//                h_Mag_Array    : Host array storing the prepared B field   data of each target patch
//                h_Pot_Array    : Host array storing the prepared potential data of each target patch
//                h_Corner_Array : Array storing the physical corner coordinates of each patch
//                NPatchGroup    : Number of patch groups evaluated simultaneously by GPU
//                dh             : Grid size
//                Safety         : dt safety factor
//                MinPres        : Minimum allowed pressure
//                P5_Gradient    : Use 5-points stencil to evaluate the potential gradient
//                UsePot         : Add self-gravity and/or external potential
//                ExtAcc         : Add external acceleration
//                TargetTime     : Target physical time
//                GPU_NStream    : Number of CUDA streams for the asynchronous memory copy
//
// Return      :  h_dt_Array
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_dtSolver( const Solver_t TSolver, real h_dt_Array[], const real h_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                          const real h_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const real h_Pot_Array[][ CUBE(GRA_NXT) ],
                          const double h_Corner_Array[][3], const int NPatchGroup, const real dh, const real Safety,
                          const real MinPres, const bool P5_Gradient, const bool UsePot, const OptExtAcc_t ExtAcc,
                          const double TargetTime, const int GPU_NStream )
{

// check
#  ifdef GAMER_DEBUG
   if ( TSolver != DT_FLU_SOLVER )
#  ifdef GRAVITY
   if ( TSolver != DT_GRA_SOLVER )
#  endif
      Aux_Error( ERROR_INFO, "TSolver != DT_FLU_SOLVER / DT_GRA_SOLVER !!\n" );

   if ( h_dt_Array == NULL )
      Aux_Error( ERROR_INFO, "h_dt_Array == NULL !!\n" );

   if ( TSolver == DT_FLU_SOLVER  &&  h_Flu_Array == NULL )
      Aux_Error( ERROR_INFO, "h_Flu_Array == NULL !!\n" );

#  ifdef GRAVITY
   if ( TSolver == DT_GRA_SOLVER )
   {
      if ( UsePot  &&  h_Pot_Array == NULL )
         Aux_Error( ERROR_INFO, "h_Pot_Array == NULL !!\n" );

      if ( ExtAcc )
      {
         if ( h_Corner_Array     == NULL )   Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
         if ( d_Corner_Array_PGT == NULL )   Aux_Error( ERROR_INFO, "d_Corner_Array_PGT == NULL !!\n" );
      }
   }
#  endif

#  ifdef MHD
   if ( TSolver == DT_FLU_SOLVER  &&  h_Mag_Array == NULL )
      Aux_Error( ERROR_INFO, "h_Mag_Array == NULL !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG


// set the block size
   const int NPatch = NPatchGroup*8;
   dim3 BlockDim_dtSolver( 1, 1, 1 );

   switch ( TSolver )
   {
      case DT_FLU_SOLVER:
         BlockDim_dtSolver.x = DT_FLU_BLOCK_SIZE;
      break;

#     ifdef GRAVITY
      case DT_GRA_SOLVER:
         BlockDim_dtSolver.x = DT_GRA_BLOCK_SIZE;
      break;
#     endif

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
   }


// set the number of patches and the corresponding data size to be transferred into GPU in each stream
   int *NPatch_per_Stream = new int [GPU_NStream];
   int *UsedPatch         = new int [GPU_NStream];
   int *dt_MemSize        = new int [GPU_NStream];
   int *Corner_MemSize    = new int [GPU_NStream];

   int *Flu_MemSize       = ( TSolver == DT_FLU_SOLVER ) ? new int [GPU_NStream] : NULL;
#  ifdef MHD
   int *Mag_MemSize       = ( TSolver == DT_FLU_SOLVER ) ? new int [GPU_NStream] : NULL;
#  endif
#  ifdef GRAVITY
   int *Pot_MemSize       = ( TSolver == DT_GRA_SOLVER ) ? new int [GPU_NStream] : NULL;
#  endif


// number of patches in each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatch;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatch / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatch - UsedPatch[GPU_NStream-1];
   }

// corresponding data size to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      switch ( TSolver )
      {
         case DT_FLU_SOLVER:
            Flu_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*CUBE(PS1)*FLU_NIN_T;
#           ifdef MHD
            Mag_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*PS1P1*SQR(PS1)*NCOMP_MAG;
#           endif
         break;

#        ifdef GRAVITY
         case DT_GRA_SOLVER:
            Pot_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*CUBE(GRA_NXT);
            Corner_MemSize[s] = sizeof(double)*NPatch_per_Stream[s]*3;
         break;
#        endif

         default :
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
      }

      dt_MemSize[s] = sizeof(real)*NPatch_per_Stream[s];
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      switch ( TSolver )
      {
         case DT_FLU_SOLVER:
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_T      + UsedPatch[s], h_Flu_Array    + UsedPatch[s],
                               Flu_MemSize[s],    cudaMemcpyHostToDevice, Stream[s] )  );
#           ifdef MHD
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Mag_Array_T      + UsedPatch[s], h_Mag_Array    + UsedPatch[s],
                               Mag_MemSize[s],    cudaMemcpyHostToDevice, Stream[s] )  );
#           endif
         break;

#        ifdef GRAVITY
         case DT_GRA_SOLVER:
            if ( UsePot )
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_T      + UsedPatch[s], h_Pot_Array    + UsedPatch[s],
                               Pot_MemSize[s],    cudaMemcpyHostToDevice, Stream[s] )  );

            if ( ExtAcc )
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_PGT + UsedPatch[s], h_Corner_Array + UsedPatch[s],
                               Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
         break;
#        endif

         default :
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
      }
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if   ( MODEL == HYDRO )
      switch ( TSolver )
      {
         case DT_FLU_SOLVER:
            CUFLU_dtSolver_HydroCFL <<< NPatch_per_Stream[s], BlockDim_dtSolver, 0, Stream[s] >>>
                                    ( d_dt_Array_T  + UsedPatch[s],
                                      d_Flu_Array_T + UsedPatch[s],
                                      d_Mag_Array_T + UsedPatch[s],
                                      dh, Safety, MinPres,
                                      EoS_DensEint2Pres_GPUPtr, EoS_DensPres2CSqr_GPUPtr, EoS_Temper2CSqr_GPUPtr,
                                      EoS_GuessHTilde_GPUPtr, EoS_HTilde2Temp_GPUPtr );
         break;

#        ifdef GRAVITY
         case DT_GRA_SOLVER:
            CUPOT_dtSolver_HydroGravity <<< NPatch_per_Stream[s], BlockDim_dtSolver, 0, Stream[s] >>>
                                        ( d_dt_Array_T       + UsedPatch[s],
                                          d_Pot_Array_T      + UsedPatch[s],
                                          d_Corner_Array_PGT + UsedPatch[s],
                                          dh, Safety, P5_Gradient, UsePot, ExtAcc, GPUExtAcc_Ptr, TargetTime );
         break;
#        endif

         default :
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
      }

#     elif ( MODEL == ELBDM )

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

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_dt_Array + UsedPatch[s], d_dt_Array_T + UsedPatch[s],
                         dt_MemSize[s], cudaMemcpyDeviceToHost, Stream[s] )  );
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] dt_MemSize;
   delete [] Corner_MemSize;
   delete [] Flu_MemSize;
#  ifdef MHD
   delete [] Mag_MemSize;
#  endif
#  ifdef GRAVITY
   delete [] Pot_MemSize;
#  endif

} // FUNCTION : CUAPI_Asyn_dtSolver



#endif // #ifdef GPU
