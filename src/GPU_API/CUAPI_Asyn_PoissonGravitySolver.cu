#include "CUAPI.h"
#include "CUPOT.h"

#if ( defined GPU  &&  defined GRAVITY )



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

#elif ( POT_SCHEME == MG  )
__global__ void CUPOT_PoissonSolver_MG( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                        const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                              real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                        const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                                        const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                                        const IntScheme_t IntScheme );
#endif // POT_SCHEME

__global__
void CUPOT_ExtPotSolver( real g_Pot_Array[][ CUBE(GRA_NXT) ],
                         const double g_Corner_Array[][3],
                         const real g_ExtPotTable[],
                         const real dh, const ExtPot_t ExtPot_Func,
                         const double Time, const bool PotIsInit );


// Gravity solver prototypes
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


// declare all device pointers
extern real (*d_Rho_Array_P    )[ CUBE(RHO_NXT) ];
extern real (*d_Pot_Array_P_In )[ CUBE(POT_NXT) ];
extern real (*d_Pot_Array_P_Out)[ CUBE(GRA_NXT) ];
extern real (*d_Flu_Array_G    )[GRA_NIN][ CUBE(PS1)];
extern double (*d_Corner_Array_PGT)[3];
#if ( MODEL == HYDRO )
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_G)[ CUBE(USG_NXT_G) ];
extern real (*d_Flu_Array_USG_G)[GRA_NIN-1][ CUBE(PS1) ];
#else
static real (*d_Pot_Array_USG_G)[ CUBE(USG_NXT_G) ] = NULL;
static real (*d_Flu_Array_USG_G)[GRA_NIN-1][ CUBE(PS1) ] = NULL;
#endif
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_G)[ CUBE(PS1) ];
#else
static char (*d_DE_Array_G)[ CUBE(PS1) ] = NULL;
#endif
#ifdef MHD
extern real (*d_Emag_Array_G)[ CUBE(PS1) ];
#else
static real (*d_Emag_Array_G)[ CUBE(PS1) ] = NULL;
#endif
#endif // #if ( MODEL == HYDRO )
extern real  *d_ExtPotTable;

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_PoissonGravitySolver
// Description :  Invoke the CUPOT_PoissonSolver_XXtoXXcube and/or CUPOT_GravitySolver kernel(s) to evaluate
//                the gravitational potential and/or advance the fluid variables by the gravitational
//                acceleration for a group of patches
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
// Parameter   :  h_Rho_Array        : Host array storing the input density
//                h_Pot_Array_In     : Host array storing the input "coarse-grid" potential for interpolation
//                h_Pot_Array_Out    : Host array to store the output potential
//                h_Flu_Array        : Host array to store the fluid variables for the Gravity solver
//                h_Corner_Array     : Host array storing the physical corner coordinates of each patch
//                h_Pot_Array_USG    : Host array storing the prepared potential for UNSPLIT_GRAVITY
//                h_Flu_Array_USG    : Host array storing the prepared density + momentum for UNSPLIT_GRAVITY
//                h_DE_Array         : Host array storing the dual-energy status (for both input and output)
//                h_Emag_Array       : Host array storing the cell-centered magnetic energy (MHD only)
//                NPatchGroup        : Number of patch groups evaluated simultaneously by GPU
//                dt                 : Time interval to advance solution
//                dh                 : Grid size
//                SOR_Min_Iter       : Minimum # of iterations for SOR
//                SOR_Max_Iter       : Maximum # of iterations for SOR
//                SOR_Omega          : Over-relaxation parameter
//                MG_Max_Iter        : Maximum number of iterations for multigrid
//                MG_NPre_Smooth     : Number of pre-smoothing steps for multigrid
//                MG_NPos_tSmooth    : Number of post-smoothing steps for multigrid
//                MG_Tolerated_Error : Maximum tolerated error for multigrid
//                Poi_Coeff          : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                IntScheme          : Interpolation scheme for potential
//                                     --> currently supported schemes include
//                                         INT_CQUAD : conservative quadratic interpolation
//                                         INT_QUAD  : quadratic interpolation
//                P5_Gradient        : Use 5-points stencil to evaluate the potential gradient
//                ELBDM_Eta          : Particle mass / Planck constant in ELBDM
//                ELBDM_Lambda       : Quartic self-interaction coefficient in ELBDM
//                Poisson            : true --> compute the self-gravity potential and/or external potential
//                GraAcc             : true --> compute the gravitational acceleration (which can include
//                                              self-gravity, external potential, and external acceleration)
//                                              to update fluid
//                SelfGravity        : Add self-gravity potential
//                ExtPot             : Add external potential
//                ExtAcc             : Add external acceleration
//                TimeNew            : Physical time at the current  step (for external gravity)
//                TimeOld            : Physical time at the previous step (for external gravity in UNSPLIT_GRAVITY)
//                MinEint            : Internal energy floor
//                GPU_NStream        : Number of CUDA streams for the asynchronous memory copy
//
// Useless parameters in HYDRO : ELBDM_Eta, ELBDM_Lambda
// Useless parameters in ELBDM : P5_Gradient
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                                      const real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                            real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                            real h_Flu_Array    [][GRA_NIN][PS1][PS1][PS1],
                                      const double h_Corner_Array[][3],
                                      const real h_Pot_Array_USG[][USG_NXT_G][USG_NXT_G][USG_NXT_G],
                                      const real h_Flu_Array_USG[][GRA_NIN-1][PS1][PS1][PS1],
                                            char h_DE_Array     [][PS1][PS1][PS1],
                                      const real h_Emag_Array   [][PS1][PS1][PS1],
                                      const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter,
                                      const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                                      const int MG_NPre_Smooth, const int MG_NPost_Smooth,
                                      const real MG_Tolerated_Error, const real Poi_Coeff,
                                      const IntScheme_t IntScheme, const bool P5_Gradient, const real ELBDM_Eta,
                                      const real ELBDM_Lambda, const bool Poisson, const bool GraAcc,
                                      const bool SelfGravity, const OptExtPot_t ExtPot, const OptExtAcc_t ExtAcc,
                                      const double TimeNew, const double TimeOld, const real MinEint,
                                      const int GPU_NStream )
{

// model-independent constants
#  if   ( POT_SCHEME == SOR )
   const dim3 Poi_Block_Dim( RHO_NXT/2, RHO_NXT, POT_BLOCK_SIZE_Z );
#  elif ( POT_SCHEME == MG )
   const dim3 Poi_Block_Dim( POT_BLOCK_SIZE_X, 1, 1 );
#  endif
   const dim3 ExtPot_Block_Dim( EXTPOT_BLOCK_SIZE );
   const dim3 Gra_Block_Dim( GRA_BLOCK_SIZE );
   const int  NPatch      = NPatchGroup*8;
#  if   ( POT_SCHEME == SOR )
   const real Poi_Const   = Poi_Coeff*dh*dh;
   const real SOR_Omega_6 = SOR_Omega/6.0;
#  endif

// model-dependent constants
#  if   ( MODEL == HYDRO )

#  elif ( MODEL == ELBDM )
   const real ELBDM_EtaDt = ELBDM_Eta*dt;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


// check
#  if ( MODEL == ELBDM  &&  !defined STORE_POT_GHOST  &&  GRA_GHOST_SIZE != 0 )
#  warning : WARNING : GRA_GHOST_SIZE != 0 in ELBDM (without STORE_POT_GHOST) !!
#  endif

#  ifdef GAMER_DEBUG
   const int Poi_NThread = Poi_Block_Dim.x * Poi_Block_Dim.y * Poi_Block_Dim.z;

   if ( Poisson )
   {
      if ( SelfGravity )
      {
//       minimum number of threads for spatial interpolation
         if ( Poi_NThread < (POT_NXT-2)*(POT_NXT-2) )
            Aux_Error( ERROR_INFO, "Poi_NThread (%d) < (POT_NXT-2)*(POT_NXT-2) (%d) !!\n",
                       Poi_NThread, (POT_NXT-2)*(POT_NXT-2) );

//       constraint due to the reduction operation in "CUPOT_Poisson_10to14cube" and "CUPOT_PoissonSolver_MG"
#        if (  ( POT_SCHEME == SOR && defined USE_PSOLVER_10TO14 )  ||  POT_SCHEME == MG  )
         if ( Poi_NThread < 64 )
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (must >= 64) !!\n", "Poi_NThread", Poi_NThread );
#        endif

//       constraint in "CUPOT_PoissonSolver_SOR_16to18cube"
#        if ( POT_SCHEME == SOR  &&  !defined USE_PSOLVER_10TO14 )
         if ( Poi_NThread != RHO_NXT*RHO_NXT/2 )
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (must == %d) !!\n", "Poi_NThread", Poi_NThread,
                       RHO_NXT*RHO_NXT/2 );
#        endif

         if ( IntScheme != INT_CQUAD  &&  IntScheme != INT_QUAD )
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
      } // if ( SelfGravity )

      if ( ExtPot )
      {
         if ( h_Corner_Array     == NULL )   Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
         if ( d_Corner_Array_PGT == NULL )   Aux_Error( ERROR_INFO, "d_Corner_Array_PGT == NULL !!\n" );
         if ( ExtPot == EXT_POT_TABLE  &&  d_ExtPotTable == NULL )
                                             Aux_Error( ERROR_INFO, "d_ExtPotTable == NULL !!\n" );
      }
   } // if ( Poisson )

   if ( GraAcc )
   {
      if ( ExtAcc )
      {
         if ( h_Corner_Array     == NULL )   Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
         if ( d_Corner_Array_PGT == NULL )   Aux_Error( ERROR_INFO, "d_Corner_Array_PGT == NULL !!\n" );
      }

#     ifdef UNSPLIT_GRAVITY
      if ( SelfGravity  ||  ExtPot )
      {
         if ( h_Pot_Array_USG   == NULL )    Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );
         if ( d_Pot_Array_USG_G == NULL )    Aux_Error( ERROR_INFO, "d_Pot_Array_USG_G == NULL !!\n" );
      }

      if ( h_Flu_Array_USG   == NULL )       Aux_Error( ERROR_INFO, "h_Flu_Array_USG == NULL !!\n" );
      if ( d_Flu_Array_USG_G == NULL )       Aux_Error( ERROR_INFO, "d_Flu_Array_USG_G == NULL !!\n" );
#     endif

#     ifdef DUAL_ENERGY
      if ( h_DE_Array   == NULL )            Aux_Error( ERROR_INFO, "h_DE_Array == NULL !!\n" );
      if ( d_DE_Array_G == NULL )            Aux_Error( ERROR_INFO, "d_DE_Array_G == NULL !!\n" );
#     endif

#     ifdef MHD
      if ( h_Emag_Array   == NULL )          Aux_Error( ERROR_INFO, "h_Emag_Array == NULL !!\n" );
      if ( d_Emag_Array_G == NULL )          Aux_Error( ERROR_INFO, "d_Emag_Array_G == NULL !!\n" );
#     endif
   } // if ( GraAcc )
#  endif // #ifdef GAMER_DEBUG


   int *NPatch_per_Stream = new int [GPU_NStream];
   int *Rho_MemSize       = new int [GPU_NStream];
   int *Pot_MemSize_In    = new int [GPU_NStream];
   int *Pot_MemSize_Out   = new int [GPU_NStream];
   int *Flu_MemSize       = new int [GPU_NStream];
   int *Corner_MemSize    = new int [GPU_NStream];
   int *UsedPatch         = new int [GPU_NStream];
#  ifdef UNSPLIT_GRAVITY
   int *Pot_USG_MemSize   = new int [GPU_NStream];
   int *Flu_USG_MemSize   = new int [GPU_NStream];
#  endif
#  ifdef DUAL_ENERGY
   int *DE_MemSize        = new int [GPU_NStream];
#  endif
#  ifdef MHD
   int *Emag_MemSize      = new int [GPU_NStream];
#  endif


// set the number of patches in each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatch;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatch/GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatch - UsedPatch[GPU_NStream-1];
   }


// set the size of data to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Rho_MemSize    [s] = NPatch_per_Stream[s]*CUBE(RHO_NXT  )*sizeof(real);
      Pot_MemSize_In [s] = NPatch_per_Stream[s]*CUBE(POT_NXT  )*sizeof(real);
      Pot_MemSize_Out[s] = NPatch_per_Stream[s]*CUBE(GRA_NXT  )*sizeof(real);
      Flu_MemSize    [s] = NPatch_per_Stream[s]*CUBE(PS1      )*sizeof(real)*GRA_NIN;
      Corner_MemSize [s] = NPatch_per_Stream[s]*3              *sizeof(double);
#     ifdef UNSPLIT_GRAVITY
      Pot_USG_MemSize[s] = NPatch_per_Stream[s]*CUBE(USG_NXT_G)*sizeof(real);
      Flu_USG_MemSize[s] = NPatch_per_Stream[s]*CUBE(PS1      )*sizeof(real)*(GRA_NIN-1);
#     endif
#     ifdef DUAL_ENERGY
      DE_MemSize     [s] = NPatch_per_Stream[s]*CUBE(PS1      )*sizeof(char);
#     endif
#     ifdef MHD
      Emag_MemSize   [s] = NPatch_per_Stream[s]*CUBE(PS1      )*sizeof(real);
#     endif
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      if ( Poisson )
      {
         if ( SelfGravity )
         {
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Rho_Array_P      + UsedPatch[s], h_Rho_Array    + UsedPatch[s],
                                                Rho_MemSize[s],    cudaMemcpyHostToDevice, Stream[s] )  );

            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_P_In   + UsedPatch[s], h_Pot_Array_In + UsedPatch[s],
                                                Pot_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
         }

         if ( ExtPot )
         {
            CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_PGT + UsedPatch[s], h_Corner_Array + UsedPatch[s],
                                                Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
         }
      } // if ( Poisson )

      if ( GraAcc )
      {
//       no need to transfer potential if we are also invoking the Poisson solver, for which
//       potential data will already be in GPU
         if (  ( SelfGravity || ExtPot )  &&  !Poisson  )
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_P_Out  + UsedPatch[s], h_Pot_Array_Out + UsedPatch[s],
                                             Pot_MemSize_Out[s], cudaMemcpyHostToDevice, Stream[s] )  );

         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_G      + UsedPatch[s], h_Flu_Array     + UsedPatch[s],
                                             Flu_MemSize[s],     cudaMemcpyHostToDevice, Stream[s] )  );

         if ( ExtAcc )
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_PGT + UsedPatch[s], h_Corner_Array  + UsedPatch[s],
                                             Corner_MemSize[s],  cudaMemcpyHostToDevice, Stream[s] )  );
#        ifdef UNSPLIT_GRAVITY
         if ( SelfGravity  ||  ExtPot )
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_USG_G  + UsedPatch[s], h_Pot_Array_USG + UsedPatch[s],
                                             Pot_USG_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );

         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_USG_G  + UsedPatch[s], h_Flu_Array_USG + UsedPatch[s],
                                             Flu_USG_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
#        endif

#        ifdef DUAL_ENERGY
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_DE_Array_G      + UsedPatch[s], h_DE_Array      + UsedPatch[s],
                                             DE_MemSize[s],      cudaMemcpyHostToDevice, Stream[s] )  );
#        endif

#        ifdef MHD
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Emag_Array_G    + UsedPatch[s], h_Emag_Array    + UsedPatch[s],
                                             Emag_MemSize[s],    cudaMemcpyHostToDevice, Stream[s] )  );
#        endif
      } // if ( GraAcc )
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

//    b1. Poisson solver
      if ( Poisson )
      {
//       b1-1. self-gravity
         if ( SelfGravity )
         {
#           if ( POT_SCHEME == SOR )

#           ifdef USE_PSOLVER_10TO14
            CUPOT_PoissonSolver_SOR_10to14cube <<< NPatch_per_Stream[s], Poi_Block_Dim, 0, Stream[s] >>>
                                               ( d_Rho_Array_P     + UsedPatch[s],
                                                 d_Pot_Array_P_In  + UsedPatch[s],
                                                 d_Pot_Array_P_Out + UsedPatch[s],
                                                 SOR_Min_Iter, SOR_Max_Iter, SOR_Omega_6, Poi_Const, IntScheme );
#           else
            CUPOT_PoissonSolver_SOR_16to18cube <<< NPatch_per_Stream[s], Poi_Block_Dim, 0, Stream[s] >>>
                                               ( d_Rho_Array_P     + UsedPatch[s],
                                                 d_Pot_Array_P_In  + UsedPatch[s],
                                                 d_Pot_Array_P_Out + UsedPatch[s],
                                                 SOR_Min_Iter, SOR_Max_Iter, SOR_Omega_6, Poi_Const, IntScheme );
#           endif // #ifdef USE_PSOLVER_10TO14 ... else ...

#           elif ( POT_SCHEME == MG  )

            CUPOT_PoissonSolver_MG             <<< NPatch_per_Stream[s], Poi_Block_Dim, 0, Stream[s] >>>
                                               ( d_Rho_Array_P     + UsedPatch[s],
                                                 d_Pot_Array_P_In  + UsedPatch[s],
                                                 d_Pot_Array_P_Out + UsedPatch[s],
                                                 dh, MG_Max_Iter, MG_NPre_Smooth, MG_NPost_Smooth, MG_Tolerated_Error,
                                                 Poi_Coeff, IntScheme );

#           else

#           error : unsupported GPU Poisson solver

#           endif // POT_SCHEME
         } // if ( SelfGravity )

//       b1-2. external potential
         if ( ExtPot )
         {
            CUPOT_ExtPotSolver <<< NPatch_per_Stream[s], ExtPot_Block_Dim, 0, Stream[s] >>>
                               ( d_Pot_Array_P_Out  + UsedPatch[s],
                                 d_Corner_Array_PGT + UsedPatch[s],
                                 d_ExtPotTable,
                                 dh, GPUExtPot_Ptr, TimeNew, SelfGravity );
         }
      } // if ( Poisson )


//    b2. Gravity solver
      if ( GraAcc )
      {
#        if   ( MODEL == HYDRO )
         CUPOT_HydroGravitySolver <<< NPatch_per_Stream[s], Gra_Block_Dim, 0, Stream[s] >>>
                                  ( d_Flu_Array_G      + UsedPatch[s],
                                    d_Pot_Array_P_Out  + UsedPatch[s],
                                    d_Corner_Array_PGT + UsedPatch[s],
                                    d_Pot_Array_USG_G  + UsedPatch[s],
                                    d_Flu_Array_USG_G  + UsedPatch[s],
                                    d_DE_Array_G       + UsedPatch[s],
                                    d_Emag_Array_G     + UsedPatch[s],
                                    dt, dh, P5_Gradient,
                                    (SelfGravity || ExtPot), ExtAcc, GPUExtAcc_Ptr,
                                    TimeNew, TimeOld, MinEint,
                                    EoS_GuessHTilde_GPUPtr,
                                    EoS_Temp2HTilde_GPUPtr,
                                    EoS_HTilde2Temp_GPUPtr );

#        elif ( MODEL == ELBDM )
         CUPOT_ELBDMGravitySolver <<< NPatch_per_Stream[s], Gra_Block_Dim, 0, Stream[s] >>>
                                  ( d_Flu_Array_G      + UsedPatch[s],
                                    d_Pot_Array_P_Out  + UsedPatch[s],
                                    d_Corner_Array_PGT + UsedPatch[s],
                                    ELBDM_EtaDt, dh, ELBDM_Lambda );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL
      } // if ( GraAcc )

      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      if ( Poisson  &&  ( SelfGravity || ExtPot )  )
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Pot_Array_Out + UsedPatch[s], d_Pot_Array_P_Out + UsedPatch[s],
                                             Pot_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( GraAcc )
      {
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array     + UsedPatch[s], d_Flu_Array_G     + UsedPatch[s],
                                             Flu_MemSize[s],     cudaMemcpyDeviceToHost, Stream[s] )  );

#        ifdef DUAL_ENERGY
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_DE_Array      + UsedPatch[s], d_DE_Array_G      + UsedPatch[s],
                                             DE_MemSize[s],      cudaMemcpyDeviceToHost, Stream[s] )  );
#        endif
      }
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] Rho_MemSize;
   delete [] Pot_MemSize_In;
   delete [] Pot_MemSize_Out;
   delete [] Flu_MemSize;
   delete [] Corner_MemSize;
   delete [] UsedPatch;
#  ifdef UNSPLIT_GRAVITY
   delete [] Pot_USG_MemSize;
   delete [] Flu_USG_MemSize;
#  endif
#  ifdef DUAL_ENERGY
   delete [] DE_MemSize;
#  endif
#  ifdef MHD
   delete [] Emag_MemSize;
#  endif

} // FUNCTION : CUAPI_Asyn_PoissonGravitySolver



#endif // #if ( defined GPU  &&  defined GRAVITY )
