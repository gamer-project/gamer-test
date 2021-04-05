#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ResetParameter
// Description :  Reset runtime parameters
//
// Note        :  1. Parameters are reset here usually because they are either non-deterministic when
//                   calling Init_Load_Parameter() (because they depend on compilation options and/or other
//                   runtime parameters) or are useless/unsupported in the adopted compilation options
//                2. This function must be invoked AFTER both Init_Load_Parameter() and Init_Unit()
//                   --> The latter may reset several physical constants (e.g., ELBDM_MASS) which are
//                       required here
//                3. This function also sets the default values for the derived runtime parameters
//                   (i.e., those depend on the input runtime parameters, e.g., amr->dh[]/BoxSize[]/BoxScale[])
//                4. This function also converts all input parameters to the code units if necessary
//                   --> Only for those not in code units when loading from the input parameter files
//                       (e.g., SF_CREATE_STAR_MIN_GAS_DENS & SF_CREATE_STAR_MIN_STAR_MASS)
//-------------------------------------------------------------------------------------------------------
void Init_ResetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// helper macro for printing warning messages
#  define FORMAT_INT    %- 21d
#  define FORMAT_FLT    %- 21.14e
#  define PRINT_WARNING( name, format, reason )                                                                \
   {                                                                                                           \
      if ( MPI_Rank == 0 )                                                                                     \
         Aux_Message( stderr, "WARNING : parameter [%-28s] is reset to [" EXPAND_AND_QUOTE(format) "] %s\n",   \
                      #name, name, reason );                                                                   \
   }


// number of OpenMP threads
#  ifdef OPENMP
   if ( OMP_NTHREAD <= 0 )
   {
      OMP_NTHREAD = omp_get_max_threads();

      PRINT_WARNING( OMP_NTHREAD, FORMAT_INT, "" );
   }
#  else
   if ( OMP_NTHREAD != 1 )
   {
      OMP_NTHREAD = 1;

      PRINT_WARNING( OMP_NTHREAD, FORMAT_INT, "since OPENMP is disabled" );
   }
#  endif // #ifdef OPENMP ... else ...


// fluid dt
   if ( DT__FLUID < 0.0 )
   {
#     if   ( MODEL == HYDRO )
#     if   ( FLU_SCHEME == RTVD )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == MHM )
      DT__FLUID = 0.80;
#     elif ( FLU_SCHEME == MHM_RP )
      DT__FLUID = 0.80;
#     elif ( FLU_SCHEME == CTU )
      DT__FLUID = 0.50;
#     else
#     error : unsupported CPU hydro scheme
#     endif

#     elif  ( MODEL == ELBDM )
#     ifdef GRAVITY
      DT__FLUID = 0.20;                   // 1D k-max mode rotates 0.20*2*PI
#     else
#     ifdef LAPLACIAN_4TH
      DT__FLUID = SQRT(27.0)*M_PI/32.0;   // stability limit (~0.51)
#     else
      DT__FLUID = SQRT(3.0)*M_PI/8.0;     // stability limit (~0.68)
#     endif
#     endif // #ifdef GRAVITY ... else ...

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      PRINT_WARNING( DT__FLUID, FORMAT_FLT, "" );
   } // if ( DT__FLUID < 0.0 )

   if ( DT__FLUID_INIT < 0.0 )
   {
      DT__FLUID_INIT = DT__FLUID;

      PRINT_WARNING( DT__FLUID_INIT, FORMAT_FLT, "" );
   }


// gravity dt
#  ifdef GRAVITY
   if ( DT__GRAVITY < 0.0 )
   {
#     if   ( MODEL == HYDRO )
      DT__GRAVITY = 0.50;
#     elif  ( MODEL == ELBDM )
      DT__GRAVITY = 0.125;
#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      PRINT_WARNING( DT__GRAVITY, FORMAT_FLT, "" );
   } // if ( DT__GRAVITY < 0.0 )
#  endif


// particle dt
#  if ( defined PARTICLE  &&  !defined STORE_PAR_ACC )
   if ( DT__PARACC != 0.0 )
   {
      DT__PARACC = 0.0;    // disable it

      PRINT_WARNING( DT__PARACC, FORMAT_FLT, "since STORE_PAR_ACC is disabled" );
   }
#  endif


// Poisson solver parameters
#  ifdef GRAVITY
#  if   ( POT_SCHEME == SOR )
   Init_Set_Default_SOR_Parameter( SOR_OMEGA, SOR_MAX_ITER, SOR_MIN_ITER );
#  elif ( POT_SCHEME == MG  )
   Init_Set_Default_MG_Parameter( MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, MG_TOLERATED_ERROR );
#  endif
#  endif // GRAVITY


// external potential table
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_TABLE  &&  EXT_POT_TABLE_FLOAT8 < 0 )
   {
//    set EXT_POT_TABLE_FLOAT8 = FLOAT8 by default
#     ifdef FLOAT8
      EXT_POT_TABLE_FLOAT8 = 1;
#     else
      EXT_POT_TABLE_FLOAT8 = 0;
#     endif

      PRINT_WARNING( EXT_POT_TABLE_FLOAT8, FORMAT_INT, "to be consistent with FLOAT8" );
   }
#  endif


// GPU parameters when using CPU only (must set OMP_NTHREAD in advance)
#  ifndef GPU
   GPU_NSTREAM = 1;

   PRINT_WARNING( GPU_NSTREAM, FORMAT_INT, "since GPU is disabled" );

   if ( FLU_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      FLU_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      FLU_GPU_NPGROUP = 1;
#     endif

      PRINT_WARNING( FLU_GPU_NPGROUP, FORMAT_INT, "since GPU is disabled" );
   }

#  ifdef GRAVITY
   if ( POT_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      POT_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      POT_GPU_NPGROUP = 1;
#     endif

      PRINT_WARNING( POT_GPU_NPGROUP, FORMAT_INT, "since GPU is disabled" );
   }
#  endif

#  ifdef SUPPORT_GRACKLE
   if ( CHE_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      CHE_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      CHE_GPU_NPGROUP = 1;
#     endif

      PRINT_WARNING( CHE_GPU_NPGROUP, FORMAT_INT, "since GPU is disabled" );
   }
#  endif

   if ( SRC_GPU_NPGROUP <= 0 )
   {
#     ifdef OPENMP
      SRC_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      SRC_GPU_NPGROUP = 1;
#     endif

      PRINT_WARNING( SRC_GPU_NPGROUP, FORMAT_INT, "since GPU is disabled" );
   }
#  endif // #ifndef GPU


// derived parameters related to the simulation scale
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)     amr->dh[lv] = BOX_SIZE / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)
   {
      amr->BoxSize  [d] = NX0_TOT[d]*amr->dh   [0];
      amr->BoxScale [d] = NX0_TOT[d]*amr->scale[0];
      amr->BoxEdgeL [d] = 0.0;
      amr->BoxEdgeR [d] = amr->BoxSize[d];
      amr->BoxCenter[d] = 0.5*( amr->BoxEdgeL[d] + amr->BoxEdgeR[d] );
   }


// workload weighting at each level
// --> treat OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE the same as DT_LEVEL_DIFF_BY_2 here since we don't know the number of updates
//     at each level during the initialization
   for (int lv=0; lv<NLEVEL; lv++)
      amr->NUpdateLv[lv] = ( OPT__DT_LEVEL == DT_LEVEL_SHARED ) ? 1L : (1L<<lv);


// whether of not to allocate fluxes at the coarse-fine boundaries
#  if   ( MODEL == HYDRO )
   if ( OPT__FIXUP_FLUX )  amr->WithFlux = true;
#  elif ( MODEL == ELBDM )
   if ( OPT__FIXUP_FLUX )  amr->WithFlux = true;
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// whether of not to allocate electric field arrays at the coarse-fine boundaries
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )    amr->WithElectric = true;
#  endif


// ELBDM parameters
#  if ( MODEL == ELBDM )
   ELBDM_ETA = ELBDM_MASS / ELBDM_PLANCK_CONST;

   PRINT_WARNING( ELBDM_ETA, FORMAT_FLT, "" );

#  ifdef COMOVING
   if ( MPI_Rank == 0 )
   {
      const double JeansK = pow( 6.0*A_INIT*SQR(ELBDM_ETA), 0.25 );
      Aux_Message( stderr, "          --> corresponding initial Jean's wavenumber (wavelength) = %13.7e h/Mpc (%13.7e Mpc/h))\n",
                   JeansK, 2.0*M_PI/JeansK );
   }
#  endif

   if ( ELBDM_TAYLOR3_AUTO )
   {
      ELBDM_TAYLOR3_COEFF = NULL_REAL;

      PRINT_WARNING( ELBDM_TAYLOR3_COEFF, FORMAT_FLT, "since ELBDM_TAYLOR3_AUTO is enabled" );
   }
#  endif // #if ( MODEL == ELBDM )


// interpolation schemes for the fluid variables
#  if   ( MODEL == HYDRO )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_CQUAD;

      PRINT_WARNING( OPT__FLU_INT_SCHEME, FORMAT_INT, "" );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_CQUAD;

      PRINT_WARNING( OPT__REF_FLU_INT_SCHEME, FORMAT_INT, "" );
   }

#  elif ( MODEL == ELBDM )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_CQUAR;

      PRINT_WARNING( OPT__FLU_INT_SCHEME, FORMAT_INT, "" );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_CQUAR;

      PRINT_WARNING( OPT__REF_FLU_INT_SCHEME, FORMAT_INT, "" );
   }

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// initial dump ID
   if ( INIT_DUMPID < 0 )  DumpID = 0;
   else                    DumpID = INIT_DUMPID;


// ResPower2 in the AMR_t structure
   int NBits0, NX0_TOT_Max;

   NX0_TOT_Max = ( NX0_TOT[0] > NX0_TOT[1]  ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_TOT_Max = ( NX0_TOT[2] > NX0_TOT_Max ) ? NX0_TOT[2] : NX0_TOT_Max;
   NBits0      = (int)log2( (double)NX0_TOT_Max );

   if (  ( NX0_TOT_Max & (NX0_TOT_Max-1) ) != 0  )    NBits0 ++;  // check if NX0_TOT_Max is a power of two

   for (int lv=0; lv<NLEVEL; lv++)  amr->ResPower2[lv] = NBits0 + lv;


// particle options
#  ifdef PARTICLE
// check if the periodic BC is applied to all directions
   bool PeriodicAllDir = true;
   for (int t=0; t<6; t++)
   {
      if ( OPT__BC_FLU[t] != BC_FLU_PERIODIC )
      {
         PeriodicAllDir = false;
         break;
      }
   }

// set RemoveCell to the distance where potential extrapolation is required when adopting non-periodic BC
   if ( !PeriodicAllDir  &&  amr->Par->RemoveCell < 0.0 )
   {
      switch ( amr->Par->Interp )
      {
         case ( PAR_INTERP_NGP ):   amr->Par->RemoveCell = 1.0;   break;
         case ( PAR_INTERP_CIC ):   amr->Par->RemoveCell = 1.5;   break;
         case ( PAR_INTERP_TSC ):   amr->Par->RemoveCell = 2.0;   break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }

      const double PAR_REMOVE_CELL = amr->Par->RemoveCell;
      PRINT_WARNING( PAR_REMOVE_CELL, FORMAT_FLT, "for the adopted PAR_INTERP scheme" );
   }

// RemoveCell is useless for the periodic B.C.
   else if ( PeriodicAllDir  &&  amr->Par->RemoveCell >= 0.0 )
   {
      amr->Par->RemoveCell = -1.0;

      const double PAR_REMOVE_CELL = amr->Par->RemoveCell;
      PRINT_WARNING( PAR_REMOVE_CELL, FORMAT_FLT, "since the periodic BC is adopted along all directions" );
   }

// number of ghost zones for the particle interpolation scheme
   if ( amr->Par->GhostSize < 0 )
   {
      switch ( amr->Par->Interp )
      {
         case ( PAR_INTERP_NGP ): amr->Par->GhostSize = 0;  break;
         case ( PAR_INTERP_CIC ): amr->Par->GhostSize = 1;  break;
         case ( PAR_INTERP_TSC ): amr->Par->GhostSize = 1;  break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }

      PRINT_WARNING( amr->Par->GhostSize, FORMAT_INT, "for the adopted PAR_INTERP scheme" );
   }
#  endif // #ifdef PARTICLE


// Green's function coefficient at the origin (must be set after setting amr->Par->Interp)
#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_ISOLATED  &&  GFUNC_COEFF0 < 0.0 )
   {
      /*
#     ifdef PARTICLE
      switch ( amr->Par->Interp )
      {
         case ( PAR_INTERP_NGP ):   GFUNC_COEFF0 = 0.0;   break;
         case ( PAR_INTERP_CIC ):   GFUNC_COEFF0 = 4.0;   break;
         case ( PAR_INTERP_TSC ):   GFUNC_COEFF0 = 4.8;   break;
         default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      }
#     else
      GFUNC_COEFF0 = 0.0;
#     endif
      */

      GFUNC_COEFF0 = 3.8;  // empirically determined value for minimizing the center-of-mass drift

      PRINT_WARNING( GFUNC_COEFF0, FORMAT_FLT, "" );
   }
#  endif


// 1st-order flux correction
#  if ( MODEL == HYDRO )
   if ( OPT__1ST_FLUX_CORR < 0 )
   {
#     ifdef MHD
      OPT__1ST_FLUX_CORR = FIRST_FLUX_CORR_NONE;

      PRINT_WARNING( OPT__1ST_FLUX_CORR, FORMAT_INT, "for MHD" );

#     else

#     if ( FLU_SCHEME == RTVD )
      OPT__1ST_FLUX_CORR = FIRST_FLUX_CORR_NONE;
#     else
      OPT__1ST_FLUX_CORR = FIRST_FLUX_CORR_3D1D;
#     endif

      PRINT_WARNING( OPT__1ST_FLUX_CORR, FORMAT_INT, "for HYDRO" );
#     endif // #ifdef MHD ... else ...
   }

   if      ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_NONE  &&  OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_NONE )
   {
      OPT__1ST_FLUX_CORR_SCHEME = RSOLVER_1ST_NONE;

      PRINT_WARNING( OPT__1ST_FLUX_CORR_SCHEME, FORMAT_INT, "since OPT__1ST_FLUX_CORR is disabled" );
   }

   else if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE  &&  OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_DEFAULT )
   {
#     ifdef MHD
      OPT__1ST_FLUX_CORR_SCHEME = RSOLVER_1ST_HLLE;
#     else
      OPT__1ST_FLUX_CORR_SCHEME = RSOLVER_1ST_HLLE;
#     endif

      PRINT_WARNING( OPT__1ST_FLUX_CORR_SCHEME, FORMAT_INT, "" );
   }
#  endif // if ( MODEL == HYDRO )


// timing options
   if ( OPT__TIMING_BARRIER < 0 )
   {
#     ifdef TIMING
      if ( OPT__TIMING_BALANCE  ||  OPT__TIMING_MPI )
         OPT__TIMING_BARRIER = 1;
      else
#     endif
         OPT__TIMING_BARRIER = 0;

#     ifdef TIMING_SOLVER
      OPT__TIMING_BARRIER = 1;
#     endif

      PRINT_WARNING( OPT__TIMING_BARRIER, FORMAT_INT, "" );
   }


// physical time
   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef COMOVING
      Time     [lv] = A_INIT;          // will be overwritten during restart
#     else
      Time     [lv] = 0.0;             // will be overwritten during restart
#     endif
      Time_Prev[lv] = -__FLT_MAX__;    // initialize as negative to indicate that it has not been set yet

      amr->FluSgTime[lv][   amr->FluSg[lv] ] = Time[lv];
#     ifdef MHD
      amr->MagSgTime[lv][   amr->MagSg[lv] ] = Time[lv];
#     endif
#     ifdef GRAVITY
      amr->PotSgTime[lv][   amr->PotSg[lv] ] = Time[lv];
#     endif

      amr->FluSgTime[lv][ 1-amr->FluSg[lv] ] = Time_Prev[lv];
#     ifdef MHD
      amr->MagSgTime[lv][ 1-amr->MagSg[lv] ] = Time_Prev[lv];
#     endif
#     ifdef GRAVITY
      amr->PotSgTime[lv][ 1-amr->PotSg[lv] ] = Time_Prev[lv];
#     endif
   }


// OPT__CORR_AFTER_ALL_SYNC
   if ( OPT__CORR_AFTER_ALL_SYNC == CORR_AFTER_SYNC_DEFAULT )
   {
      OPT__CORR_AFTER_ALL_SYNC = CORR_AFTER_SYNC_BEFORE_DUMP;

      PRINT_WARNING( OPT__CORR_AFTER_ALL_SYNC, FORMAT_INT, "" );
   }


// turn off "OPT__OVERLAP_MPI" if (1) OVERLAP_MPI=ff, (2) SERIAL=on, (3) LOAD_BALANCE=off,
//                                (4) OPENMP=off, (5) MPI thread support=MPI_THREAD_SINGLE
#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      PRINT_WARNING( OPT__OVERLAP_MPI, FORMAT_INT, "since OVERLAP_MPI is disabled in the makefile" );
   }
#  endif

#  ifdef SERIAL
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      PRINT_WARNING( OPT__OVERLAP_MPI, FORMAT_INT, "since SERIAL is enabled" );
   }
#  endif // ifdef SERIAL

#  ifndef LOAD_BALANCE
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      PRINT_WARNING( OPT__OVERLAP_MPI, FORMAT_INT, "since LOAD_BALANCE is disabled" );
   }
#  endif // #ifndef LOAD_BALANCE

#  ifndef OPENMP
   if ( OPT__OVERLAP_MPI )
   {
      OPT__OVERLAP_MPI = false;

      PRINT_WARNING( OPT__OVERLAP_MPI, FORMAT_INT, "since OPENMP is disabled" );
   }
#  endif

#  ifndef SERIAL
// check the level of MPI thread support
   int MPI_Thread_Status;
   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
   {
      OPT__OVERLAP_MPI = false;

      PRINT_WARNING( OPT__OVERLAP_MPI, FORMAT_INT, "since the level of MPI thread support == MPI_THREAD_SINGLE" );
   }
#  endif


// disable "OPT__CK_FLUX_ALLOCATE" if no flux arrays are going to be allocated
   if ( OPT__CK_FLUX_ALLOCATE  &&  !amr->WithFlux )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      PRINT_WARNING( OPT__CK_FLUX_ALLOCATE, FORMAT_INT, "since no flux is required" );
   }


// no temporal interpolation in the shared time-step integration
   if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  OPT__INT_TIME )
   {
      OPT__INT_TIME = false;

      PRINT_WARNING( OPT__INT_TIME, FORMAT_INT, "since OPT__DT_LEVEL == DT_LEVEL_SHARED" );
   }


// always turn on "OPT__VERBOSE" in the debug mode
#  ifdef GAMER_DEBUG
   if ( !OPT__VERBOSE )
   {
      OPT__VERBOSE = true;

      PRINT_WARNING( OPT__VERBOSE, FORMAT_INT, "since GAMER_DEBUG is enabled" );
   }
#  endif


// flux operations are useful in HYDRO/ELBDM only
#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
   if ( OPT__FIXUP_FLUX )
   {
      OPT__FIXUP_FLUX = false;

      PRINT_WARNING( OPT__FIXUP_FLUX, FORMAT_INT, "since it's only supported in HYDRO/ELBDM" );
   }

   if ( OPT__CK_FLUX_ALLOCATE )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      PRINT_WARNING( OPT__CK_FLUX_ALLOCATE, FORMAT_INT, "since it's only supported in HYDRO/ELBDM" );
   }
#  endif


// turn off refinement criteria and checks related to density if "DENS" is not defined
#  ifndef DENS
   if ( OPT__FLAG_RHO )
   {
      OPT__FLAG_RHO = false;

      PRINT_WARNING( OPT__FLAG_RHO, FORMAT_INT, "since the symbolic constant DENS is not defined" );
   }

   if ( OPT__FLAG_RHO_GRADIENT )
   {
      OPT__FLAG_RHO_GRADIENT = false;

      PRINT_WARNING( OPT__FLAG_RHO_GRADIENT, FORMAT_INT, "since the symbolic constant DENS is not defined" );
   }

   if ( OPT__CK_REFINE )
   {
      OPT__CK_REFINE = false;

      PRINT_WARNING( OPT__CK_REFINE, FORMAT_INT, "since the symbolic constant DENS is not defined" );
   }
#  endif // #ifndef DENS


// conservation check is supported only in HYDRO/ELBDM
#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
   if ( OPT__CK_CONSERVATION )
   {
      OPT__CK_CONSERVATION = false;

      PRINT_WARNING( OPT__CK_CONSERVATION, FORMAT_INT, "since it's only supported in HYDRO/ELBDM" );
   }
#  endif


// disable OPT__LR_LIMITER if it is useless
#  if ( MODEL == HYDRO  &&  FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
   if ( OPT__LR_LIMITER != LR_LIMITER_NONE )
   {
      OPT__LR_LIMITER = LR_LIMITER_NONE;

      PRINT_WARNING( OPT__LR_LIMITER, FORMAT_INT, "since it's only useful for the MHM/MHM_RP/CTU schemes" );
   }
#  endif


// disable the refinement flag of Jeans length if GRAVITY is disabled
#  if ( MODEL == HYDRO  &&  !defined GRAVITY )
   if ( OPT__FLAG_JEANS )
   {
      OPT__FLAG_JEANS = false;

      PRINT_WARNING( OPT__FLAG_JEANS, FORMAT_INT, "since GRAVITY is disabled" );
   }
#  endif


// flux operation in ELBDM is useful only if CONSERVE_MASS is on
#  if ( MODEL == ELBDM  &&  !defined CONSERVE_MASS )
   if ( OPT__FIXUP_FLUX )
   {
      OPT__FIXUP_FLUX = false;

      PRINT_WARNING( OPT__FIXUP_FLUX, FORMAT_INT, "since CONSERVE_MASS is disabled" );
   }

   if ( OPT__CK_FLUX_ALLOCATE )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      PRINT_WARNING( OPT__CK_FLUX_ALLOCATE, FORMAT_INT, "since CONSERVE_MASS is disabled" );
   }
#  endif


// OPT__OUTPUT_BASEPS is not supported if GRAVITY is disabled
#  ifndef GRAVITY
   if ( OPT__OUTPUT_BASEPS )
   {
      OPT__OUTPUT_BASEPS = false;

      PRINT_WARNING( OPT__OUTPUT_BASEPS, FORMAT_INT, "since GRAVITY is disabled" );
   }
#  endif


// reset MPI_NRank_X
#  ifdef SERIAL
   for (int d=0; d<3; d++)
   if ( MPI_NRank_X[d] != 1 )
   {
      MPI_NRank_X[d] = 1;
      PRINT_WARNING( MPI_NRank_X[d], FORMAT_INT, "for SERIAL" );
   }
#  endif

#  ifdef LOAD_BALANCE
   for (int d=0; d<3; d++)
   if ( MPI_NRank_X[d] > 0 )
   {
      MPI_NRank_X[d] = -1;

      PRINT_WARNING( MPI_NRank_X[d], FORMAT_INT, "since it's useless" );
   }
#  endif


// turn off various timing options if TIMING is disabled
#  ifndef TIMING
   if ( OPT__RECORD_PERFORMANCE )
   {
      OPT__RECORD_PERFORMANCE = false;

      PRINT_WARNING( OPT__RECORD_PERFORMANCE, FORMAT_INT, "since TIMING is disabled" );
   }

   if ( OPT__TIMING_BARRIER != 0 )
   {
      OPT__TIMING_BARRIER = 0;

      PRINT_WARNING( OPT__TIMING_BARRIER, FORMAT_INT, "since TIMING is disabled" );
   }

   if ( OPT__TIMING_BALANCE )
   {
      OPT__TIMING_BALANCE = false;

      PRINT_WARNING( OPT__TIMING_BALANCE, FORMAT_INT, "since TIMING is disabled" );
   }

   if ( OPT__TIMING_MPI )
   {
      OPT__TIMING_MPI = false;

      PRINT_WARNING( OPT__TIMING_MPI, FORMAT_INT, "since TIMING is disabled" );
   }
#  endif // #ifndef TIMING


// only load-balance routines support OPT__TIMING_MPI
#  ifndef LOAD_BALANCE
   if ( OPT__TIMING_MPI )
   {
      OPT__TIMING_MPI = false;

      PRINT_WARNING( OPT__TIMING_MPI, FORMAT_INT, "since LOAD_BALANCE is disabled" );
   }
#  endif


// OPT__UM_IC_NVAR
   if ( OPT__INIT == INIT_BY_FILE  &&  OPT__UM_IC_NVAR <= 0 )
   {
#     if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
      OPT__UM_IC_NVAR = NCOMP_TOTAL - 1;  // do not load the dual-energy field from the disk

#     elif ( MODEL == ELBDM )
      OPT__UM_IC_NVAR = NCOMP_TOTAL - 1;  // do not load the density field from the disk

#     else
      OPT__UM_IC_NVAR = NCOMP_TOTAL;      // load all fields
#     endif

      PRINT_WARNING( OPT__UM_IC_NVAR, FORMAT_INT, "" );
   }


// always turn on "OPT__CK_PARTICLE" when debugging particles
#  ifdef DEBUG_PARTICLE
   if ( !OPT__CK_PARTICLE )
   {
      OPT__CK_PARTICLE = true;

      PRINT_WARNING( OPT__CK_PARTICLE, FORMAT_INT, "since DEBUG_PARTICLE is enabled" );
   }
#  endif


// set particle initialization to PAR_INIT_BY_RESTART for restart
#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_RESTART  &&  amr->Par->Init != PAR_INIT_BY_RESTART )
   {
      amr->Par->Init = PAR_INIT_BY_RESTART;

      const ParInit_t PAR_INIT = amr->Par->Init;
      PRINT_WARNING( PAR_INIT, FORMAT_INT, "for restart" );
   }
#  endif


// JEANS_MIN_PRES must work with GRAVITY
#  if ( MODEL == HYDRO )
#  ifndef GRAVITY
   if ( JEANS_MIN_PRES )
   {
      JEANS_MIN_PRES = false;

      PRINT_WARNING( JEANS_MIN_PRES, FORMAT_INT, "since GRAVITY is disabled" );
   }
#  endif

   if ( JEANS_MIN_PRES  &&  JEANS_MIN_PRES_LEVEL < 0 )
   {
      JEANS_MIN_PRES_LEVEL = MAX_LEVEL;

      PRINT_WARNING( JEANS_MIN_PRES_LEVEL, FORMAT_INT, "" );
   }
#  endif


// MIN_PRES and MIN_EINT
#  if ( MODEL == HYDRO )
   if      ( MIN_PRES > 0.0  &&  MIN_EINT == 0.0 )
   {
      MIN_EINT = MIN_PRES*1.5;

      PRINT_WARNING( MIN_EINT, FORMAT_FLT, "" );
   }

   else if ( MIN_EINT > 0.0  &&  MIN_PRES == 0.0 )
   {
      MIN_PRES = MIN_EINT/1.5;

      PRINT_WARNING( MIN_PRES, FORMAT_FLT, "" );
   }
#  endif


// AUTO_REDUCE_DT only works for DT_LEVEL_FLEXIBLE
   if ( AUTO_REDUCE_DT  &&  OPT__DT_LEVEL != DT_LEVEL_FLEXIBLE )
   {
      AUTO_REDUCE_DT = false;

      PRINT_WARNING( AUTO_REDUCE_DT, FORMAT_INT, "since OPT__DT_LEVEL != DT_LEVEL_FLEXIBLE" );
   }


// FLAG_BUFFER_SIZE at the level MAX_LEVEL-1 and MAX_LEVEL-2
   if ( FLAG_BUFFER_SIZE_MAXM1_LV < 0 )
   {
      FLAG_BUFFER_SIZE_MAXM1_LV = FLAG_BUFFER_SIZE;

      PRINT_WARNING( FLAG_BUFFER_SIZE_MAXM1_LV, FORMAT_INT, "" );
   }

// must set FLAG_BUFFER_SIZE_MAXM1_LV in advance
   if ( FLAG_BUFFER_SIZE_MAXM2_LV < 0 )
   {
      FLAG_BUFFER_SIZE_MAXM2_LV = ( FLAG_BUFFER_SIZE_MAXM1_LV + FLAG_BUFFER_SIZE_MAXM1_LV%2 + PS1 ) / 2;

      PRINT_WARNING( FLAG_BUFFER_SIZE_MAXM2_LV, FORMAT_INT, "" );
   }


// star-formation options
#  ifdef STAR_FORMATION
   if ( SF_CREATE_STAR_MIN_LEVEL < 0 )
   {
      SF_CREATE_STAR_MIN_LEVEL = MAX_LEVEL;

      PRINT_WARNING( SF_CREATE_STAR_MIN_LEVEL, FORMAT_INT, "" );
   }

   if ( SF_CREATE_STAR_DET_RANDOM < 0 )
   {
#     ifdef BITWISE_REPRODUCIBILITY
         SF_CREATE_STAR_DET_RANDOM = 1;
         PRINT_WARNING( SF_CREATE_STAR_DET_RANDOM, FORMAT_INT, "since BITWISE_REPRODUCIBILITY is enabled" );
#     else
         SF_CREATE_STAR_DET_RANDOM = 0;
         PRINT_WARNING( SF_CREATE_STAR_DET_RANDOM, FORMAT_INT, "since BITWISE_REPRODUCIBILITY is disabled" );
#     endif

   }
#  endif // #ifdef STAR_FORMATION


// feedback options
#  ifdef FEEDBACK
   if ( FB_LEVEL < 0 )
   {
      FB_LEVEL = MAX_LEVEL;

      PRINT_WARNING( FB_LEVEL, FORMAT_INT, "" );
   }
#  endif // #ifdef FEEDBACK


// convert to code units
#  ifdef STAR_FORMATION
// SF_CREATE_STAR_MIN_GAS_DENS: HI count/cm^3 --> mass density in code units
   SF_CREATE_STAR_MIN_GAS_DENS *= Const_mH / UNIT_D;

   PRINT_WARNING( SF_CREATE_STAR_MIN_GAS_DENS, FORMAT_FLT, "to be consistent with the code units" );


// SF_CREATE_STAR_MIN_STAR_MASS: Msun --> code units
   SF_CREATE_STAR_MIN_STAR_MASS *= Const_Msun / UNIT_M;

   PRINT_WARNING( SF_CREATE_STAR_MIN_STAR_MASS, FORMAT_FLT, "to be consistent with the code units" );
#  endif // #ifdef STAR_FORMATION


// disable OPT__MINIMIZE_MPI_BARRIER in the serial mode
#  ifdef SERIAL
   if ( OPT__MINIMIZE_MPI_BARRIER )
   {
      OPT__MINIMIZE_MPI_BARRIER = false;

      PRINT_WARNING( OPT__MINIMIZE_MPI_BARRIER, FORMAT_INT, "since SERIAL is enabled" );
   }
#  endif


// disable OPT__INIT_GRID_WITH_OMP if OPENMP is disabled
#  ifndef OPENMP
   if ( OPT__INIT_GRID_WITH_OMP )
   {
      OPT__INIT_GRID_WITH_OMP = false;

      PRINT_WARNING( OPT__INIT_GRID_WITH_OMP, FORMAT_INT, "since OPENMP is disabled" );
   }
#  endif


// remove symbolic constants and macros only used in this structure
#  undef FORMAT_INT
#  undef FORMAT_FLT
#  undef QUOTE


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ResetParameter
