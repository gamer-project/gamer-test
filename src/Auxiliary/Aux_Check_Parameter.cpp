#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Parameter
// Description :  Check the initial parameter setting
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Parameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ...\n" );


// general errors
// =======================================================================================
#  if ( PATCH_SIZE%2 != 0 )
#     error : ERROR : PATCH_SIZE must be an even number !!
#  endif

#  if ( defined TIMING_SOLVER  &&  !defined TIMING )
#     error : ERROR : TIMING_SOLVER must work with TIMING !!
#  endif

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP; the macro "_OPENMP" is NOT defined !!
#  endif

#  if ( defined OVERLAP_MPI  &&  !defined LOAD_BALANCE )
#     error : ERROR : OVERLAP_MPI must work with LOAD_BALANCE !!
#  endif

#  if ( !defined GRAVITY  &&  defined UNSPLIT_GRAVITY )
#     error : ERROR : UNSPLIT_GRAVITY must work with GRAVITY !!
#  endif

#  if ( defined UNSPLIT_GRAVITY  &&  MODEL != HYDRO )
#     error : ERROR : currently UNSPLIT_GRAVITY is only supported in HYDRO !!
#  endif

#  if ( NCOMP_PASSIVE < 0 )
#     error : ERROR : incorrect number of NCOMP_PASSIVE !!
#  endif

#  ifdef SERIAL
   int NRank = 1;
#  else
   int NRank, MPI_Thread_Status, MPI_Init_Status;

   MPI_Initialized( &MPI_Init_Status );
   if ( MPI_Init_Status == false )  Aux_Error( ERROR_INFO, "MPI_Init() has not been called !!\n" );

   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT supported since the level of MPI thread support == %s\n",
                 "OPT__OVERLAP_MPI", "MPI_THREAD_SINGLE" );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
#  endif

   if ( MPI_NRank != NRank )
      Aux_Error( ERROR_INFO, "MPI_NRank (%d) != MPI_Comm_size (%d) !!\n", MPI_NRank, NRank );

   for (int d=0; d<3; d++)
      if ( NX0_TOT[d]%PS2 != 0 )
         Aux_Error( ERROR_INFO, "NX0_TOT_%c (%d) is NOT a multiple of %d (i.e., two patches) !!\n", 'X'+d, NX0_TOT[d], PS2 );

   if ( END_STEP < 0  &&  OPT__INIT != INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d\" [>=0] !!\n", "END_STEP", END_STEP );

   if ( END_T < 0.0  &&  OPT__INIT != INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %14.7e\" [>=0] !!\n", "END_T", END_T );

#  ifndef LOAD_BALANCE
   if ( NX0_TOT[0]%(PS2*MPI_NRank_X[0]) != 0  ||  NX0_TOT[1]%(PS2*MPI_NRank_X[1]) != 0  ||
        NX0_TOT[2]%(PS2*MPI_NRank_X[2]) != 0 )
      Aux_Error( ERROR_INFO, "number of base-level patches in each direction and in each MPI rank must be \"%s\" !!\n",
                 "a multiple of TWO" );

   if ( MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank )
      Aux_Error( ERROR_INFO, "MPI_NRANK_X (%d) * MPI_NRANK_Y (%d) * MPI_NRANK_Z (%d) = %d != MPI_Comm_size (%d) !!\n",
                 MPI_NRank_X[0], MPI_NRank_X[1], MPI_NRank_X[2], MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2], MPI_NRank );
#  endif

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_STEP  &&  OUTPUT_STEP <= 0 )
      Aux_Error( ERROR_INFO, "OUTPUT_STEP (%ld) <= 0 !!\n", OUTPUT_STEP );

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT  &&  OUTPUT_DT <= 0.0 )
      Aux_Error( ERROR_INFO, "OUTPUT_DT (%14.7e) <= 0.0 !!\n", OUTPUT_DT );

#  ifndef SUPPORT_HDF5
   if ( OPT__OUTPUT_TOTAL == OUTPUT_FORMAT_HDF5 )
      Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for OPT__OUTPUT_TOTAL == 1 !!\n" );
#  endif

   if (  ( OPT__OUTPUT_PART == OUTPUT_YZ  ||  OPT__OUTPUT_PART == OUTPUT_Y  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_X < 0.0  ||  OUTPUT_PART_X >= amr->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_X (out of range [0<=X<%lf]) !!\n", amr->BoxSize[0] );

   if (  ( OPT__OUTPUT_PART == OUTPUT_XZ  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_Y < 0.0  ||  OUTPUT_PART_Y >= amr->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Y (out of range [0<=Y<%lf]) !!\n", amr->BoxSize[1] );

   if (  ( OPT__OUTPUT_PART == OUTPUT_XY  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Y )  &&
         ( OUTPUT_PART_Z < 0.0  ||  OUTPUT_PART_Z >= amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Z (out of range [0<=Z<%lf]) !!\n", amr->BoxSize[2] );

   if (  OPT__OUTPUT_PART == OUTPUT_DIAG  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "\"%s\" only works with CUBIC domain !!\n",
                 "OPT__OUTPUT_PART == 7 (OUTPUT_DIAG)" );

   if (  OPT__OUTPUT_BASEPS  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "\"%s\" only works with CUBIC domain !!\n", "OPT__OUTPUT_BASEPS" );

   if ( OPT__CK_REFINE  &&  !OPT__FLAG_RHO )
      Aux_Error( ERROR_INFO, "currently the check \"%s\" must work with \"%s\" !!\n",
                 "OPT__CK_REFINE", "OPT__FLAG_RHO" );

#  if   ( MODEL == HYDRO )
   if (  ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP )
         &&  Flu_ParaBuf < 2  )
      Aux_Error( ERROR_INFO, "Lohner error estimator does NOT work when Flu_ParaBuf (%d) < 2 !!\n", Flu_ParaBuf );
#  elif ( MODEL == ELBDM )
   if (  OPT__FLAG_LOHNER_DENS  &&  Flu_ParaBuf < 2  )
      Aux_Error( ERROR_INFO, "Lohner error estimator does NOT work when Flu_ParaBuf (%d) < 2 !!\n", Flu_ParaBuf );
#  else
#  error : ERROR : unsupported MODEL !!
#  endif

#  ifdef GPU
#  ifdef LAOHU
   if ( OPT__GPUID_SELECT < -3 )
#  else
   if ( OPT__GPUID_SELECT < -2 )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n", "OPT__GPUID_SELECT", OPT__GPUID_SELECT );
#  endif
#  endif

#  ifdef SERIAL
   if ( MPI_NRank != 1 )   Aux_Error( ERROR_INFO, "\"MPI_NRank != 1\" in the serial code !!\n" );

   for (int d=0; d<3; d++)
      if ( MPI_NRank_X[d] != 1 )
         Aux_Error( ERROR_INFO, "\"MPI_NRank_%c (%d) != 1\" in the serial code !!\n", 'X'+d, MPI_NRank_X[d] );
#  endif // #ifdef SERIAL

#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT enabled in the Makefile for \"%s\" !!\n",
                 "OVERLAP_MPI", "OPT__OVERLAP_MPI" );
#  endif

   if ( OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT supported yet !!\n", "OPT__OVERLAP_MPI" );

   if ( AUTO_REDUCE_DT )
   {
      if ( OPT__OVERLAP_MPI )
         Aux_Error( ERROR_INFO, "currently \"%s\" does not work with \"%s\" !!\n",
                    "AUTO_REDUCE_DT", "OPT__OVERLAP_MPI" );

      if ( OPT__DT_LEVEL != DT_LEVEL_FLEXIBLE )
         Aux_Error( ERROR_INFO, "\"%s\" must work with \"%s\" !!\n",
                    "AUTO_REDUCE_DT", "OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE" );
   }

#  if ( MODEL != HYDRO )
   for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_OUTFLOW )
         Aux_Error( ERROR_INFO, "outflow boundary condition (OPT__BC_FLU=2) only works with HYDRO !!\n" );

   for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_REFLECTING )
         Aux_Error( ERROR_INFO, "reflecting boundary condition (OPT__BC_FLU=3) only works with HYDRO !!\n" );
#  endif

   for (int f=0; f<6; f+=2)
      if (  ( OPT__BC_FLU[f] == BC_FLU_PERIODIC  &&  OPT__BC_FLU[f+1] != BC_FLU_PERIODIC )  ||
            ( OPT__BC_FLU[f] != BC_FLU_PERIODIC  &&  OPT__BC_FLU[f+1] == BC_FLU_PERIODIC )   )
         Aux_Error( ERROR_INFO, "periodic and non-periodic boundary conditions cannot be mixed along the same dimension"
                                "--> please modify OPT__BC_FLU[%d/%d] !!\n", f, f+1 );

#  ifndef TIMING
   if ( OPT__TIMING_MPI )  Aux_Error( ERROR_INFO, "OPT__TIMING_MPI must work with TIMING !!\n" );
#  endif

   if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  OPT__INT_TIME )
      Aux_Error( ERROR_INFO, "OPT__INT_TIME should be disabled when \"OPT__DT_LEVEL == DT_LEVEL_SHARED\" !!\n" );

   if ( INT_MONO_COEFF < 1.0  ||  INT_MONO_COEFF > 4.0 )
      Aux_Error( ERROR_INFO, "INT_MONO_COEFF (%14.7e) is not within the correct range [1.0, 4.0] !!\n", INT_MONO_COEFF );

   if ( OPT__MEMORY_POOL  &&  !OPT__REUSE_MEMORY )
      Aux_Error( ERROR_INFO, "please turn on OPT__REUSE_MEMORY for OPT__MEMORY_POOL !!\n" );

   if ( OPT__CORR_AFTER_ALL_SYNC != CORR_AFTER_SYNC_NONE  &&  OPT__CORR_AFTER_ALL_SYNC != CORR_AFTER_SYNC_EVERY_STEP  &&
        OPT__CORR_AFTER_ALL_SYNC != CORR_AFTER_SYNC_BEFORE_DUMP )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__CORR_AFTER_ALL_SYNC = %d\" [0/1/2] !!\n", OPT__CORR_AFTER_ALL_SYNC );

   if ( OPT__MINIMIZE_MPI_BARRIER )
   {
#     if ( defined GRAVITY  &&  !defined STORE_POT_GHOST )
         Aux_Error( ERROR_INFO, "OPT__MINIMIZE_MPI_BARRIER must work with STORE_POT_GHOST when GRAVITY is on !!\n" );
#     endif

#     ifdef PARTICLE
      if ( !amr->Par->ImproveAcc )
         Aux_Error( ERROR_INFO, "OPT__MINIMIZE_MPI_BARRIER must work with PAR_IMPROVE_ACC when PARTICLE is on !!\n" );

#     ifdef LOAD_BALANCE
      if ( LB_INPUT__PAR_WEIGHT == 0.0  &&  MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : consider fine tuning \"LB_INPUT__PAR_WEIGHT\" to optimize the performance of OPT__MINIMIZE_MPI_BARRIER !!\n" );
#     endif
#     endif

      if ( OPT__TIMING_BARRIER )
         Aux_Error( ERROR_INFO, "OPT__MINIMIZE_MPI_BARRIER does NOT work with OPT__TIMING_BARRIER !!\n" );

      if ( AUTO_REDUCE_DT  &&  MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : AUTO_REDUCE_DT will introduce an extra MPI barrier for OPT__MINIMIZE_MPI_BARRIER !!\n" );
   }

#  ifdef BITWISE_REPRODUCIBILITY
   if ( OPT__CORR_AFTER_ALL_SYNC != CORR_AFTER_SYNC_BEFORE_DUMP  &&  OPT__CORR_AFTER_ALL_SYNC != CORR_AFTER_SYNC_EVERY_STEP )
      Aux_Error( ERROR_INFO, "please set OPT__CORR_AFTER_ALL_SYNC to 1/2 when BITWISE_REPRODUCIBILITY is enabled !!\n" );

   if ( ! OPT__FIXUP_RESTRICT )
      Aux_Error( ERROR_INFO, "must enable OPT__FIXUP_RESTRICT for BITWISE_REPRODUCIBILITY !!\n" );
#  endif

#  if ( !defined SERIAL  &&  !defined LOAD_BALANCE )
   if ( OPT__INIT == INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "must enable either SERIAL or LOAD_BALANCE for OPT__INIT=3 !!\n" );
#  endif



// general warnings
// =======================================================================================
#  ifdef OVERLAP_MPI
#     warning : WARNING : make sure to link with multi-thread supported MPI and FFTW for "OVERLAP_MPI"
#  endif

#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {
      if ( OMP_NTHREAD != omp_get_num_threads() )
         Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) != omp_get_num_threads (%d) at MPI_Rank %d !!\n",
                      OMP_NTHREAD, omp_get_num_threads(), MPI_Rank );
   }

   const int OMP_Max_NThread = omp_get_max_threads();
   if ( OMP_NTHREAD > OMP_Max_NThread )
   {
      Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) > omp_get_max_threads (%d) at MPI_Rank %d !!\n",
                   OMP_NTHREAD, OMP_Max_NThread, MPI_Rank );
   }
#  endif


   if ( MPI_Rank == 0 ) {

#  if ( defined GAMER_DEBUG  &&  !defined BITWISE_REPRODUCIBILITY )
      Aux_Message( stderr, "WARNING : you might want to turn on BITWISE_REPRODUCIBILITY for GAMER_DEBUG !!\n" );
#  endif

   if ( !OPT__OUTPUT_TOTAL  &&  !OPT__OUTPUT_PART  &&  !OPT__OUTPUT_USER  &&  !OPT__OUTPUT_BASEPS )
#  ifdef PARTICLE
   if ( !OPT__OUTPUT_PAR_TEXT )
#  endif
      Aux_Message( stderr, "WARNING : all output options are turned off --> no data will be output !!\n" );

   if ( OPT__CK_REFINE )
      Aux_Message( stderr, "WARNING : \"%s\" check may fail due to the proper-nesting constraint !!\n",
                   "OPT__CK_REFINE" );

   if ( !OPT__INIT_RESTRICT )
      Aux_Message( stderr, "WARNING : OPT__INIT_RESTRICT is disabled !!\n" );

#  ifdef TIMING_SOLVER
   Aux_Message( stderr, "WARNING : \"TIMING_SOLVER\" will disable the concurrent execution\n" );
   Aux_Message( stderr, "          between GPU and CPU and hence will decrease the overall performance !!\n" );
#  endif

   if ( OPT__CK_REFINE )
      Aux_Message( stderr, "WARNING : currently the check \"%s\" only works with \"%s\" !!\n",
                   "OPT__CK_REFINE", "OPT__FLAG_RHO == 1" );

   bool Flag = ( OPT__FLAG_RHO  ||  OPT__FLAG_RHO_GRADIENT  ||  OPT__FLAG_LOHNER_DENS  ||  OPT__FLAG_USER );
#  if ( MODEL == HYDRO )
   Flag |= OPT__FLAG_PRES_GRADIENT;
   Flag |= OPT__FLAG_ENGY_GRADIENT;
   Flag |= OPT__FLAG_VORTICITY;
   Flag |= OPT__FLAG_JEANS;
   Flag |= OPT__FLAG_LOHNER_ENGY;
   Flag |= OPT__FLAG_LOHNER_PRES;
   Flag |= OPT__FLAG_LOHNER_TEMP;
#  ifdef MHD
   Flag |= OPT__FLAG_CURRENT;
#  endif
#  ifdef SRHD
   Flag |= OPT__FLAG_4VELOCITY;
   Flag |= OPT__FLAG_MOM_OVER_DENS;
   Flag |= OPT__FLAG_LORENTZ_GRADIENT;
#  endif
#  endif
#  if ( MODEL == ELBDM )
   Flag |= OPT__FLAG_ENGY_DENSITY;
#  endif
#  ifdef PARTICLE
   Flag |= OPT__FLAG_NPAR_PATCH;
   Flag |= OPT__FLAG_NPAR_CELL;
   Flag |= OPT__FLAG_PAR_MASS_CELL;
#  endif

   if ( !Flag )
      Aux_Message( stderr, "WARNING : all flag criteria are turned off --> no refinement will be performed !!" );

   if ( OPT__NO_FLAG_NEAR_BOUNDARY  )
      Aux_Message( stderr, "WARNING : OPT__NO_FLAG_NEAR_BOUNDARY is on --> patches adjacent to the "
                           "simulation boundaries are NOT allowed for refinement !!\n" );

   if ( OPT__OVERLAP_MPI )
   {
      Aux_Message( stderr, "WARNING : \"%s\" is still experimental and is not fully optimized !!\n",
                   "OPT__OVERLAP_MPI" );

#     ifdef OPENMP
      omp_set_nested( true );

      if ( !omp_get_nested() )
         Aux_Message( stderr, "WARNING : OpenMP nested parallelism is NOT supported for \"%s\" !!\n",
                      "OPT__OVERLAP_MPI" );

      omp_set_nested( false );
#     else
      Aux_Message( stderr, "WARNING : OpenMP is NOT turned on for \"%s\" !!\n", "OPT__OVERLAP_MPI" );
#     endif
   } // if ( OPT__OVERLAP_MPI )

   if ( OPT__TIMING_BARRIER )
      Aux_Message( stderr, "WARNING : \"%s\" may deteriorate performance (especially if %s is on) ...\n",
                   "OPT__TIMING_BARRIER", "OPT__OVERLAP_MPI" );

   if ( OPT__TIMING_BARRIER  &&  !OPT__TIMING_BALANCE )
   {
      Aux_Message( stderr, "REMINDER : \"%s\" is on, but the time waiting for other ranks will NOT be included in individual timers ...\n",
                   "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "           --> the sum of individual timer may be less than the total elapsed time due to load imbalance ...\n" );
   }

#  ifdef TIMING
   if ( !OPT__TIMING_BARRIER )
   {
      Aux_Message( stderr, "REMINDER : \"%s\" is off for TIMING\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "           --> Some timing results (especially MPI and particle routines) may be less accurate due to load imbalance ...\n" );
   }
#  endif

   if (  ( OPT__TIMING_BALANCE || OPT__TIMING_MPI )  &&  !OPT__TIMING_BARRIER  )
   {
      Aux_Message( stderr, "REMINDER : \"%s\" is off for OPT__TIMING_BALANCE/OPT__TIMING_MPI\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "           --> Some timing results (especially MPI and particle routines) may be less accurate due to load imbalance ...\n" );
   }

#  ifdef PARTICLE
   if ( OPT__TIMING_BALANCE )
   {
      Aux_Message( stderr, "REMINDER : \"%s\" does NOT work well for particle routines\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "           --> Because many particle routines call MPI_Barrier implicitly\n" );
      Aux_Message( stderr, "               (so as Gra_AdvanceDt when PARTICLE is on)\n" );
   }
#  endif

#  if ( defined GRAVITY  &&  GRA_GHOST_SIZE == 0  &&  defined STORE_POT_GHOST )
   Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless when GRA_GHOST_SIZE == 0 !!\n" );
#  endif

#  if ( !defined GRAVITY  &&  defined STORE_POT_GHOST )
   Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless when GRAVITY is off !!\n" );
#  endif

#  if ( NCOMP_PASSIVE > 0 )
   if ( OPT__NORMALIZE_PASSIVE )
      Aux_Message( stderr, "REMINDER : OPT__NORMALIZE_PASSIVE will break the strict conservation of individual passive scalar\n" );
   else
   {
      Aux_Message( stderr, "REMINDER : disabling OPT__NORMALIZE_PASSIVE will break the strict equality between\n" );
      Aux_Message( stderr, "           sum(passive_scalar_mass_density) and gas_mass_density\n" );
   }
#  endif

#  if   ( MODEL == HYDRO )
   if ( ! INT_OPP_SIGN_0TH_ORDER )
      Aux_Message( stderr, "WARNING : disabling INT_OPP_SIGN_0TH_ORDER may cause unphysically large velocity during interpolation !!\n" );
#  elif ( MODEL == ELBDM )
   if (   INT_OPP_SIGN_0TH_ORDER )
      Aux_Message( stderr, "WARNING : INT_OPP_SIGN_0TH_ORDER is not recommended for ELBDM !!\n" );
#  endif

   } // if ( MPI_Rank == 0 )



// load balance
// =======================================================================================
#ifdef LOAD_BALANCE

// errors
// ------------------------------
#  ifdef SERIAL
#     error : ERROR : LOAD_BALANCE and SERIAL must NOT be enabled at the same time !!
#  endif

#  if ( LOAD_BALANCE != HILBERT )
#     error : ERROR : currently GAMER only supports "LOAD_BALANCE == HILBERT" !!
#  endif

// for sending fluid data fixed by coarse-fine fluxes correctly
   if ( OPT__FIXUP_FLUX  &&  Flu_ParaBuf >= PATCH_SIZE )
      Aux_Error( ERROR_INFO, "\"%s\" is required for \"%s\" in LOAD_BALANCE --> check LB_RecordExchangeFixUpDataPatchID() !!\n",
                 "Flu_ParaBuf < PATCH_SIZE", "OPT__FIXUP_FLUX" );

// ensure that the variable "PaddedCr1D" will not overflow
   const int Padded              = 1<<NLEVEL;
   const int BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[2]/PATCH_SIZE + 2*Padded };
   if (  (double)BoxNScale_Padded[0]*(double)BoxNScale_Padded[1]*(double)BoxNScale_Padded[2]-1.0 > (double)__ULONG_MAX__  )
   {
      Aux_Message( stderr, "ERROR : PaddedCr1D can overflow !!\n" );
      Aux_Message( stderr, "    --> The currently maximum allowed 1-D resolution is about 2^24 (for PS1=8)\n" );
      Aux_Message( stderr, "    --> Please either set NLEVEL to a smaller number or disable LOAD_BALANCE !!\n" );
      MPI_Exit();
   }


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( NX0_TOT[0] != NX0_TOT[1]  ||  NX0_TOT[0] != NX0_TOT[2] )
   {
      Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-cubic simulation box\n" );
      Aux_Message( stderr, "          (NX0_TOT[0] != NX0_TOT[1] and/or NX0_TOT[0] != NX0_TOT[2]) !!\n" );
   }

   for (int d=0; d<3; d++)
   {
      if ( NX0_TOT[d] & (NX0_TOT[d]-1) )
      {
         Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-power-of-two " );
         Aux_Message( stderr, "simulation box (NX0_TOT[%d] = %d) !!\n", d, NX0_TOT[d] );
      }
   }

   } // if ( MPI_Rank == 0 )

#else // #ifdef LOAD_BALANCE ... else ...

   if ( OPT__OVERLAP_MPI  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : currently \"%s\" is useful only in LOAD_BALANCE !!\n",
                   "OPT__OVERLAP_MPI" );

#endif // #ifdef LOAD_BALANCE ... else ...



// comoving frame
// =======================================================================================
#ifdef COMOVING

// errors
// ------------------------------
#  if ( MODEL == HYDRO )
   if ( fabs(GAMMA-5.0/3.0) > 1.0e-4 )
      Aux_Error( ERROR_INFO, "GAMMA must be equal to 5.0/3.0 for COMOVING !!\n" );
#  endif


// warnings
// ------------------------------
#  ifndef GRAVITY
   if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : \"%s\" is useless when \"%s\" is disabled !!\n",
                   "COMOVING", "GRAVITY" );
#  endif

#endif // COMOVING



// fluid solver in all models
// =======================================================================================

// errors
// ------------------------------
   if ( Flu_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Flu_ParaBuf, PATCH_SIZE );

   if ( GPU_NSTREAM < 1 )  Aux_Error( ERROR_INFO, "GPU_NSTREAM (%d) < 1 !!\n", GPU_NSTREAM );

   if ( FLU_GPU_NPGROUP % GPU_NSTREAM != 0 )
      Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 FLU_GPU_NPGROUP, GPU_NSTREAM );

#  ifdef OPENMP
   if ( FLU_GPU_NPGROUP < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", FLU_GPU_NPGROUP, OMP_NTHREAD );
#  endif

   if ( OPT__FIXUP_FLUX  &&  !amr->WithFlux )
      Aux_Error( ERROR_INFO, "%s is enabled but amr->WithFlux is off !!\n", "OPT__FIXUP_FLUX" );

#  if ( NLEVEL > 1 )
   int Trash_RefFlu, NGhost_RefFlu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, Trash_RefFlu, NGhost_RefFlu );
   if ( Flu_ParaBuf < NGhost_RefFlu )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) < NGhost_RefFlu (%d) --> refinement will fail !!\n",
                 Flu_ParaBuf, NGhost_RefFlu );
#  endif

   if ( OPT__RESET_FLUID  &&   OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT supported for \"%s\" !!\n", "OPT__OVERLAP_MPI", "OPT__RESET_FLUID" );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( OPT__CK_FLUX_ALLOCATE  &&  !amr->WithFlux )
      Aux_Message( stderr, "WARNING : \"%s\" is useless since no flux is required !!\n",
                   "OPT__CK_FLUX_ALLOCATE" );

   if ( DT__FLUID < 0.0  ||  DT__FLUID > 1.0 )
      Aux_Message( stderr, "WARNING : DT__FLUID (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__FLUID );

   if ( DT__FLUID_INIT < 0.0  ||  DT__FLUID_INIT > 1.0 )
      Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__FLUID_INIT );

   if ( OPT__RESET_FLUID  &&   OPT__INIT == INIT_BY_FILE )
      Aux_Message( stderr, "WARNING : \"%s\" will NOT be applied to the input uniform data !!\n", "OPT__RESET_FLUID" );

   } // if ( MPI_Rank == 0 )



// fluid solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( NCOMP_FLUID != 5 )
#     error : ERROR : NCOMP_FLUID != 5 in HYDRO !!
#  endif

#  if ( NCOMP_TOTAL != NFLUX_TOTAL )
#     error : ERROR : NCOMP_TOTAL != NFLUX_TOTAL !!
#  endif

#  if ( NCOMP_PASSIVE != 0  &&  FLU_SCHEME == RTVD )
#     error : RTVD does NOT support passive scalars !!
#  endif

#  if ( FLU_SCHEME != RTVD  &&  FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#     error : ERROR : unsupported hydro scheme in the makefile !!
#  endif

#  if ( defined UNSPLIT_GRAVITY  &&  FLU_SCHEME == RTVD )
#     error : ERROR : RTVD does not support UNSPLIT_GRAVITY !!
#  endif

#  if ( defined LR_SCHEME  &&  LR_SCHEME != PLM  &&  LR_SCHEME != PPM )
#     error : ERROR : unsupported data reconstruction scheme (PLM/PPM) !!
#  endif

#  ifdef MHD
#   if ( defined RSOLVER  &&  RSOLVER != ROE  &&  RSOLVER != HLLE  &&  RSOLVER != HLLD )
#     error : ERROR : unsupported Riemann solver for MHD (ROE/HLLE/HLLD) !!
#   endif
#  else
#   if ( defined RSOLVER  &&  RSOLVER != EXACT  &&  RSOLVER != ROE  &&  RSOLVER != HLLE  &&  RSOLVER != HLLC )
#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!
#   endif
#  endif // MHD

#  ifdef SRHD
#   if ( defined RSOLVER  &&  RSOLVER != HLLC  &&  RSOLVER != HLLE )
#     error : ERROR : unsupported Riemann solver for MHD (HLLC/HLLE) !!
#   endif

#  if ( defined PARTICLE )
#     error : ERROR : SRHD does not support PARTICLE !!
#  endif

#  if ( NCOMP_PASSIVE != 0 )
#     error : SRHD does NOT support passive scalars !!
#  endif

#  if ( defined DUAL_ENERGY )
#     error : ERROR : SRHD does not support DUAL_ENERGY !!
#  endif

#  if ( defined SUPPORT_GRACKLE )
#     error : ERROR : SRHD does not support GRACKLE !!
#  endif

#  if ( FLU_SCHEME != MHM && FLU_SCHEME != MHM_RP )
#     error : ERROR : unsupported FLU_SCHEME for SRHD (MHM/MHM_RP) !!
#  endif

#  if ( EOS != EOS_TAUBMATHEWS )
#     error : ERROR : unsupported EOS for SRHD (EOS_TAUBMATHEWS) !!
#  endif

#  else
#  if ( EOS == EOS_TAUBMATHEWS )
#     error : ERROR : unsupported EOS for non-SRHD (EOS_TAUBMATHEWS) !!
#  endif

#  endif // #  ifdef SRHD

#  ifdef DUAL_ENERGY
#   if ( FLU_SCHEME == RTVD )
#     error : RTVD does NOT support DUAL_ENERGY !!
#   endif

#   if ( DUAL_ENERGY != DE_ENPY )
#     error : ERROR : unsupported dual-energy formalism (DE_ENPY only, DE_EINT is not supported yet) !!
#   endif

#  if ( DUAL_ENERGY == DE_ENPY  &&  EOS != EOS_GAMMA )
#     error : ERROR : DUAL_ENERGY=DE_ENPY only supports EOS_GAMMA !!
#  endif

#   if ( DUAL_ENERGY == DE_ENPY  &&  defined COSMIC_RAY )
#     error : COSMIC_RAY does NOT support DUAL_ENERGY=DE_ENPY !!
#   endif
#  endif // #ifdef DUAL_ENERGY

#  ifdef MHD
#   if ( defined CHECK_INTERMEDIATE  &&  CHECK_INTERMEDIATE != HLLE  &&  CHECK_INTERMEDIATE != HLLD )
#     error : ERROR : unsupported option in CHECK_INTERMEDIATE (HLLE/HLLD) !!
#   endif
#  else
#   if ( defined CHECK_INTERMEDIATE  &&  CHECK_INTERMEDIATE != EXACT  &&  CHECK_INTERMEDIATE != HLLE  &&  \
        CHECK_INTERMEDIATE != HLLC )
#     error : ERROR : unsupported option in CHECK_INTERMEDIATE (EXACT/HLLE/HLLC) !!
#   endif
#  endif // MHD

#  ifdef COSMIC_RAY
#     error : ERROR : COSMIC_RAY is NOT supported yet !!
#  endif

#  if ( defined LR_EINT  &&  FLU_SCHEME == CTU )
#     error : ERROR : CTU does NOT support LR_EINT in CUFLU.h !!
#  endif

#  if ( EOS != EOS_GAMMA       &&  EOS != EOS_ISOTHERMAL  &&  EOS != EOS_NUCLEAR  && \
        EOS != EOS_TAUBMATHEWS &&  EOS != EOS_TABULAR     &&  EOS != EOS_USER )
#     error : ERROR : unsupported equation of state (EOS_GAMMA/EOS_ISOTHERMAL/EOS_NUCLEAR/EOS_TAUBMATHEWS/EOS_TABULAR/EOS_USER) !!
#  endif

#  if ( EOS != EOS_GAMMA )
#     if ( HLLC_WAVESPEED == HLL_WAVESPEED_ROE  ||  HLLE_WAVESPEED == HLL_WAVESPEED_ROE )
#        error : ERROR : HLL_WAVESPEED_ROE only works with EOS_GAMMA !!
#     endif

#     if (  defined RSOLVER  &&  ( RSOLVER == ROE || RSOLVER == EXACT )  )
#        error : ERROR : unsupported Riemann solver for EOS != EOS_GAMMA (HLLE/HLLC/HLLD) !!
#     endif

#     if ( defined LR_SCHEME  &&  defined CHAR_RECONSTRUCTION )
#        error : ERROR : CHAR_RECONSTRUCTION only works with EOS_GAMMA !!
#     endif

#     if ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU )
#        error : RTVD and CTU only support EOS_GAMMA !!
#     endif

#     ifdef COMOVING
#        error : ERROR : COMOVING currently only supports EOS_GAMMA !!
#     endif

      if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE  &&  OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_ROE )
         Aux_Error( ERROR_INFO, "OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_ROE only supports EOS_GAMMA !!\n" );

      if ( JEANS_MIN_PRES )
         Aux_Error( ERROR_INFO, "JEANS_MIN_PRES currently only supports EOS_GAMMA !!\n" );
#  endif // if ( EOS != EOS_GAMMA )

#  if ( EOS == EOS_NUCLEAR )
      Aux_Error( ERROR_INFO, "EOS_NUCLEAR is not supported yet !!\n" );
#  endif

#  if ( EOS == EOS_TABULAR )
      Aux_Error( ERROR_INFO, "EOS_TABULAR is not supported yet !!\n" );
#  endif

#  ifdef BAROTROPIC_EOS
#     if ( EOS == EOS_GAMMA  ||  EOS == EOS_NUCLEAR )
#        error : ERROR : BAROTROPIC_EOS is incompatible with EOS_GAMMA/EOS_NUCLEAR !!
#     endif
#  else
#     if ( EOS == EOS_ISOTHERMAL )
#        error : ERROR : must enable BAROTROPIC_EOS for EOS_ISOTHERMAL !!
#     endif
#  endif // #ifdef BAROTROPIC_EOS ... else ...


   if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
   {
#     ifdef MHD
      /*
      if ( OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_ROE  &&  OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_HLLD  &&
           OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_HLLE )
         Aux_Error( ERROR_INFO, "unsupported parameter \"%s = %d\" !!\n", "OPT__1ST_FLUX_CORR_SCHEME", OPT__1ST_FLUX_CORR_SCHEME );
         */
      Aux_Error( ERROR_INFO, "MHD does not support \"OPT__1ST_FLUX_CORR\" !!\n" );
#     else
      if ( OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_ROE  &&  OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_HLLC  &&
           OPT__1ST_FLUX_CORR_SCHEME != RSOLVER_1ST_HLLE )
         Aux_Error( ERROR_INFO, "unsupported parameter \"%s = %d\" !!\n", "OPT__1ST_FLUX_CORR_SCHEME", OPT__1ST_FLUX_CORR_SCHEME );
#     endif

#     if ( FLU_SCHEME == RTVD )
         Aux_Error( ERROR_INFO, "RTVD does not support \"OPT__1ST_FLUX_CORR\" !!\n" );
#     endif
   }

#  if ( FLU_SCHEME == RTVD )
   if ( JEANS_MIN_PRES )
      Aux_Error( ERROR_INFO, "RTVD does not support \"JEANS_MIN_PRES\" !!\n" );
#  endif


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( defined RSOLVER  &&  RSOLVER == EXACT )
#     warning : WARNING : exact Riemann solver is not recommended since the vacuum solution has not been implemented
      Aux_Message( stderr, "WARNING : exact Riemann solver is not recommended since the vacuum solution " );
      Aux_Message( stderr,           "has not been implemented !!\n" );
#  endif

#  if ( defined CHAR_RECONSTRUCTION  &&  defined GRAVITY )
#     warning : WARNING : "CHAR_RECONSTRUCTION" is less robust and can cause negative density/pressure !!
      Aux_Message( stderr, "WARNING : \"CHAR_RECONSTRUCTION\" is less robust and could cause negative " );
      Aux_Message( stderr,           "density/pressure !!\n" );
#  endif

   if ( !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_FLUX" );

   if ( !OPT__FIXUP_RESTRICT )
      Aux_Message( stderr, "WARNING : \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_RESTRICT" );

   if ( OPT__CK_FLUX_ALLOCATE  &&  !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : %s is useless since %s is off !!\n", "OPT__CK_FLUX_ALLOCATE", "OPT__FIXUP_FLUX" );

   if ( OPT__INIT == INIT_BY_FILE )
      Aux_Message( stderr, "WARNING : currently we don't check MIN_DENS/PRES for the initial data loaded from UM_IC !!\n" );

   if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
      Aux_Message( stderr, "REMINDER : OPT__1ST_FLUX_CORR may break the strict conservation of fluid variables\n" );

#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_ACTIVATE && OPT__FLAG_LOHNER_TEMP )
      Aux_Message( stderr, "WARNING : currently we do not use Grackle to calculate temperature for OPT__FLAG_LOHNER_TEMP !!\n" );
#  endif

   if ( ! OPT__LAST_RESORT_FLOOR )
      Aux_Message( stderr, "WARNING : disabling OPT__LAST_RESORT_FLOOR could be dangerous and is mainly for debugging only !!\n" );

   if ( MIN_DENS == 0.0 )
      Aux_Message( stderr, "WARNING : MIN_DENS == 0.0 could be dangerous and is mainly for debugging only !!\n" );
   else
      Aux_Message( stderr, "WARNING : MIN_DENS (%13.7e) is on --> please ensure that this value is reasonable !!\n", MIN_DENS );

   if ( MIN_PRES == 0.0 )
      Aux_Message( stderr, "WARNING : MIN_PRES == 0.0 could be dangerous and is mainly for debugging only !!\n" );
   else
      Aux_Message( stderr, "WARNING : MIN_PRES (%13.7e) is on --> please ensure that this value is reasonable !!\n", MIN_PRES );

   if ( MIN_EINT == 0.0 )
      Aux_Message( stderr, "WARNING : MIN_EINT == 0.0 could be dangerous and is mainly for debugging only !!\n" );
   else
      Aux_Message( stderr, "WARNING : MIN_EINT (%13.7e) is on --> please ensure that this value is reasonable !!\n", MIN_EINT );

#  if (  defined LR_EINT  &&  ( EOS == EOS_GAMMA || EOS == EOS_ISOTHERMAL )  )
      Aux_Message( stderr, "WARNING : LR_EINT is not recommended for EOS_GAMMA/EOS_ISOTHERMAL !!\n" );
#  endif
   } // if ( MPI_Rank == 0 )


// check for MHM/MHM_RP/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

// errors
// ------------------------------
#  if ( LR_SCHEME == PPM )
   if ( OPT__LR_LIMITER == EXTPRE )
      Aux_Error( ERROR_INFO, "currently the PPM reconstruction does not support the \"%s\" limiter\n",
                 "extrema-preserving" );
#  endif

   if ( OPT__LR_LIMITER != VANLEER  &&  OPT__LR_LIMITER != GMINMOD  &&  OPT__LR_LIMITER != ALBADA  &&
        OPT__LR_LIMITER != EXTPRE   &&  OPT__LR_LIMITER != VL_GMINMOD )
      Aux_Error( ERROR_INFO, "unsupported data reconstruction limiter (OPT__LR_IMITER = %d) !!\n",
                 OPT__LR_LIMITER );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   } // if ( MPI_Rank == 0 )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )


// check for MHM/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 2 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 2  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )


// check for MHM_RP
// ------------------------------
#  if ( FLU_SCHEME == MHM_RP )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM_RP )


// check for RTVD
// ------------------------------
#  if ( FLU_SCHEME == RTVD )

#  if ( FLU_GHOST_SIZE != 3 )
#     error : ERROR : please set FLU_GHOST_SIZE = 3 for the relaxing TVD scheme !!
#  endif

#  endif // if ( FLU_SCHEME == RTVD )


// check for MHD
// ------------------------------
#  ifdef MHD

#  if ( !defined SERIAL  &&  !defined LOAD_BALANCE )
#     error : ERROR : MHD must work with either SERIAL or LOAD_BALANCE !!
#  endif

#  if ( FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#     error : ERROR : unsupported MHD scheme in the makefile (MHM_RP/CTU) !!
#  endif

#  if ( HLLE_WAVESPEED == HLL_WAVESPEED_PVRS )
#     error : ERROR : HLL_WAVESPEED_PVRS does not support MHD !!
#  endif

#  if ( HLLD_WAVESPEED != HLL_WAVESPEED_DAVIS )
#     error : ERROR : HLLD_WAVESPEED only supports HLL_WAVESPEED_DAVIS !!
#  endif

   if ( OPT__MAG_INT_SCHEME != INT_MINMOD1D  &&  OPT__MAG_INT_SCHEME != INT_VANLEER  &&
        OPT__MAG_INT_SCHEME != INT_CQUAD  &&  OPT__MAG_INT_SCHEME != INT_CQUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" (2,3,4,6 only) !!\n",
                 "OPT__MAG_INT_SCHEME", OPT__MAG_INT_SCHEME );

   if ( OPT__REF_MAG_INT_SCHEME != INT_MINMOD1D  &&  OPT__REF_MAG_INT_SCHEME != INT_VANLEER  &&
        OPT__REF_MAG_INT_SCHEME != INT_CQUAD  &&  OPT__REF_MAG_INT_SCHEME != INT_CQUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" (2,3,4,6 only) !!\n",
                 "OPT__REF_MAG_INT_SCHEME", OPT__REF_MAG_INT_SCHEME );

   if ( OPT__FIXUP_ELECTRIC  &&  !amr->WithElectric )
      Aux_Error( ERROR_INFO, "%s is enabled but amr->Electric is off !!\n", "OPT__FIXUP_ELECTRIC" );

   if ( OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"OPT__OVERLAP_MPI\" is NOT supported for MHD !!\n" );

   if ( OPT__INIT == INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "MHD does NOT currently support \"OPT__INIT=3\" !!\n" );

   if ( !OPT__FIXUP_RESTRICT )
      Aux_Message( stderr, "WARNING : disabling \"OPT__FIXUP_RESTRICT\" in MHD will break the divergence-free constraint !!\n" );

   if ( !OPT__FIXUP_ELECTRIC )
      Aux_Message( stderr, "WARNING : disabling \"OPT__FIXUP_ELECTRIC\" in MHD will break the divergence-free constraint !!\n" );

   if ( !OPT__OUTPUT_CC_MAG )
      Aux_Message( stderr, "WARNING : yt requires \"OPT__OUTPUT_CC_MAG\" for analyzing magnetic field !!\n" );

#  endif // #ifdef MHD



// fluid solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )

// errors
// ------------------------------
#  if ( NCOMP_FLUID != 3 )
#     error : ERROR : NCOMP_FLUID != 3 in ELBDM !!
#  endif

#  if ( FLU_NIN != 2 )
#     error : ERROR : FLU_NIN != 2 in ELBDM !!
#  endif

#  if ( FLU_NOUT != 3 )
#     error : ERROR : FLU_NOUT != 3 in ELBDM !!
#  endif

#  if ( NCOMP_PASSIVE > 0 )
#     error : ERROR : NCOMP_PASSIVE > 0 in ELBDM (currently this model does not support passive scalars) !!
#  endif

#  ifdef QUARTIC_SELF_INTERACTION
#  ifndef GRAVITY
#     error : ERROR : currently QUARTIC_SELF_INTERACTION must work with GRAVITY !!
#  endif

#  ifdef COMOVING
#     error : ERROR : QUARTIC_SELF_INTERACTION does not work with COMOVING yet !!
#  endif
#  endif // ifdef QUARTIC_SELF_INTERACTION

   if ( ELBDM_PLANCK_CONST <= 0.0 )
      Aux_Error( ERROR_INFO, "%s (%14.7e) <= 0.0 !!\n", "ELBDM_PLANCK_CONST", ELBDM_PLANCK_CONST );

   if ( ELBDM_ETA <= 0.0 )
      Aux_Error( ERROR_INFO, "%s (%14.7e) <= 0.0 !!\n", "ELBDM_ETA", ELBDM_ETA );

   if ( OPT__INT_PHASE  &&  OPT__FLU_INT_SCHEME == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" when OPT__INT_PHASE is on !!\n",
                 "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );

   if ( MIN_DENS == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MIN_DENS == 0.0 could be dangerous and is mainly for debugging only !!\n" );
   else if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MIN_DENS (%13.7e) is on --> please ensure that this value is reasonable !!\n", MIN_DENS );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( !ELBDM_TAYLOR3_AUTO  &&  ELBDM_TAYLOR3_COEFF < 1.0/8.0 )
      Aux_Message( stderr, "WARNING : ELBDM_TAYLOR3_COEFF (%13.7e) < 0.125 is unconditionally unstable !!\n",
                   ELBDM_TAYLOR3_COEFF );

#  ifdef LAPLACIAN_4TH
   const double dt_fluid_max = 3.0*M_PI/16.0;
#  else
   const double dt_fluid_max = 0.25*M_PI;
#  endif
   if ( DT__FLUID > dt_fluid_max )
      Aux_Message( stderr, "WARNING : DT__FLUID (%13.7e) > %13.7e is unconditionally unstable (even with %s) !!\n",
                   DT__FLUID, dt_fluid_max, "ELBDM_TAYLOR3_AUTO" );

   if ( DT__FLUID_INIT > dt_fluid_max )
      Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%13.7e) > %13.7e is unconditionally unstable (even with %s) !!\n",
                   DT__FLUID_INIT, dt_fluid_max, "ELBDM_TAYLOR3_AUTO" );

   if ( !ELBDM_TAYLOR3_AUTO )
   {
//    stability limit for ELBDM_TAYLOR3_COEFF == 1.0/6.0
#     ifdef LAPLACIAN_4TH
      const double dt_fluid_max_normal = SQRT(27.0)*M_PI/32.0;
#     else
      const double dt_fluid_max_normal = SQRT(3.0)*M_PI/8.0;
#     endif

      if ( DT__FLUID > dt_fluid_max_normal  &&  ELBDM_TAYLOR3_COEFF <= 1.0/6.0 )
      {
         Aux_Message( stderr, "WARNING : DT__FLUID (%13.7e) > stability limit (%13.7e) for ELBDM_TAYLOR3_COEFF <= 1/6\n",
                      DT__FLUID, dt_fluid_max_normal );
         Aux_Message( stderr, "          --> Please either (a) set ELBDM_TAYLOR3_COEFF (%13.7e) > 1/6\n",
                      ELBDM_TAYLOR3_COEFF );
         Aux_Message( stderr, "                            (b) set DT__FLUID smaller (c) turn on ELBDM_TAYLOR3_AUTO\n" );
      }

      if ( DT__FLUID_INIT > dt_fluid_max_normal  &&  ELBDM_TAYLOR3_COEFF <= 1.0/6.0 )
      {
         Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%13.7e) > stability limit (%13.7e) for ELBDM_TAYLOR3_COEFF <= 1/6\n",
                      DT__FLUID_INIT, dt_fluid_max_normal );
         Aux_Message( stderr, "          --> Please either (a) set ELBDM_TAYLOR3_COEFF (%13.7e) > 1/6\n",
                      ELBDM_TAYLOR3_COEFF );
         Aux_Message( stderr, "                            (b) set DT__FLUID_INIT smaller (c) turn on ELBDM_TAYLOR3_AUTO\n" );
      }
   }

   if ( DT__PHASE > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PHASE (%13.7e) is not within the normal range [0...1] !!\n", DT__PHASE );

   if ( OPT__CK_FLUX_ALLOCATE  &&  !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : %s is useless since %s is off !!\n",
                   "OPT__CK_FLUX_ALLOCATE", "OPT__FIXUP_FLUX" );

#  ifdef CONSERVE_MASS
   if ( !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : %s is disabled in ELBDM even though CONSERVE_MASS is on !!\n",
                   "OPT__FIXUP_FLUX" );
#  else
   if ( OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : %s is useless in ELBDM when CONSERVE_MASS is off !!\n", "OPT__FIXUP_FLUX" );
#  endif

   if ( OPT__INIT == INIT_BY_FILE )
      Aux_Message( stderr, "WARNING : currently we don't check MIN_DENS for the initial data loaded from UM_IC !!\n" );

   } // if ( MPI_Rank == 0 )

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL



// Poisson and Gravity solvers
// =======================================================================================
#ifdef GRAVITY

// errors
// ------------------------------
#  if ( POT_SCHEME != SOR  &&  POT_SCHEME != MG )
#     error : ERROR : unsupported Poisson solver in the makefile (SOR/MG) !!
#  endif

#  if ( POT_GHOST_SIZE <= GRA_GHOST_SIZE )
      #error : ERROR : POT_GHOST_SIZE <= GRA_GHOST_SIZE !!
#  endif

#  if ( POT_GHOST_SIZE < 1 )
#     error : ERROR : POT_GHOST_SIZE < 1 !!
#  endif

#  ifdef GPU
#     if ( POT_GHOST_SIZE > 5 )
#        error : ERROR : POT_GHOST_SIZE must <= 5 for the GPU Poisson solver !!
#     endif

   if ( OPT__SELF_GRAVITY  &&  PATCH_SIZE != 8 )
      Aux_Error( ERROR_INFO, "PATCH_SIZE must == 8 for the GPU Poisson solver (OPT__SELF_GRAVITY) !!\n" );
#  endif // GPU

#  ifndef LOAD_BALANCE
   if ( NX0_TOT[2]%MPI_NRank != 0 )
   {
      Aux_Message( stderr, "ERROR : NX0_TOT[2] %% MPI_NRank != 0 !!\n" );
      Aux_Message( stderr, "--> All MPI processes must have the same number of cells in the z direction for " );
      Aux_Message( stderr,     "the slab decomposition in FFTW 2.1.5 if LOAD_BALANCE is off !!\n" );
      MPI_Exit();
   }
#  endif

   if ( Pot_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Pot_ParaBuf, PATCH_SIZE );

   if ( Rho_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Rho_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Rho_ParaBuf, PATCH_SIZE );

#  if ( POT_SCHEME == SOR )
   if ( SOR_OMEGA < 0.0 )     Aux_Error( ERROR_INFO, "SOR_OMEGA (%14.7e) < 0.0 !!\n", SOR_OMEGA );
   if ( SOR_MAX_ITER < 0 )    Aux_Error( ERROR_INFO, "SOR_MAX_ITER (%d) < 0 !!\n", SOR_MAX_ITER );
   if ( SOR_MIN_ITER < 3 )    Aux_Error( ERROR_INFO, "SOR_MIN_ITER (%d) < 3 !!\n", SOR_MIN_ITER );
#  endif

#  if ( POT_SCHEME == MG )
   if ( MG_MAX_ITER < 0 )              Aux_Error( ERROR_INFO, "MG_MAX_ITER (%d) < 0 !!\n", MG_MAX_ITER );
   if ( MG_NPRE_SMOOTH < 0 )           Aux_Error( ERROR_INFO, "MG_NPRE_SMOOTH (%d) < 0 !!\n", MG_NPRE_SMOOTH );
   if ( MG_NPOST_SMOOTH < 0 )          Aux_Error( ERROR_INFO, "MG_NPOST_SMOOTH (%d) < 0 !!\n", MG_NPOST_SMOOTH );
   if ( MG_TOLERATED_ERROR < 0.0 )     Aux_Error( ERROR_INFO, "MG_TOLERATED_ERROR (%14.7e) < 0.0 !!\n", MG_TOLERATED_ERROR );
#  endif

   if ( POT_GPU_NPGROUP % GPU_NSTREAM != 0 )
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 POT_GPU_NPGROUP, GPU_NSTREAM );

#  ifdef OPENMP
   if ( POT_GPU_NPGROUP < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", POT_GPU_NPGROUP, OMP_NTHREAD );
#  endif

#  if ( NLEVEL > 1 )
   int Trash_RefPot, NGhost_RefPot;
   Int_Table( OPT__REF_POT_INT_SCHEME, Trash_RefPot, NGhost_RefPot );
   if ( Pot_ParaBuf < NGhost_RefPot )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) < NGhost_RefPot (%d) --> refinement will fail !!\n",
                 Pot_ParaBuf, NGhost_RefPot );
#  endif

   if ( OPT__BC_POT != BC_POT_PERIODIC  &&  OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__BC_POT = %d\" [1/2] !!\n", OPT__BC_POT );

   if ( NEWTON_G <= 0.0 )     Aux_Error( ERROR_INFO, "NEWTON_G (%14.7e) <= 0.0 !!\n", NEWTON_G );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( POT_SCHEME == MG  &&  PATCH_SIZE <= 8 )
   {
      Aux_Message( stderr, "WARNING : multigrid scheme gives lower performance than SOR for " );
      Aux_Message( stderr, "PATCH_SIZE <= 8 and hence is not recommended !!\n" );
   }
#  endif

   if ( DT__GRAVITY < 0.0  ||  DT__GRAVITY > 1.0 )
      Aux_Message( stderr, "WARNING : DT__GRAVITY (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__GRAVITY );

   if ( !OPT__SELF_GRAVITY  &&  !OPT__EXT_ACC  &&  !OPT__EXT_POT )
      Aux_Message( stderr, "WARNING : all gravity options are disabled (OPT__SELF_GRAVITY, OPT__EXT_ACC, OPT__EXT_POT) !!\n" );

   } // if ( MPI_Rank == 0 )



// gravity solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( GRA_GHOST_SIZE < 1 )
#     error : ERROR : GRA_GHOST_SIZE must >= 1
#  endif

   if ( OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 1 )
      Aux_Error( ERROR_INFO, "\"%s\" requires \"%s\" !!\n",
                 "OPT__GRA_P5_GRADIENT", "GRA_GHOST_SIZE == 2" );

#  ifdef UNSPLIT_GRAVITY
   if ( OPT__GRA_P5_GRADIENT &&  USG_GHOST_SIZE_G == 1 )
      Aux_Error( ERROR_INFO, "\"%s\" requires \"%s\" for UNSPLIT_GRAVITY !!\n",
                 "OPT__GRA_P5_GRADIENT", "USG_GHOST_SIZE_G == 2" );
#  endif


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  ifndef STORE_POT_GHOST
   if ( !OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 2 )
   {
      Aux_Message( stderr, "WARNING : \"%s\" is useless when \"%s\" is disabled !!\n",
                   "GRA_GHOST_SIZE == 2", "OPT__GRA_P5_GRADIENT" );
   }
#  endif

   if ( GRA_GHOST_SIZE > 2 )  Aux_Message( stderr, "WARNING : \"GRA_GHOST_SIZE > 2\" !?\n" );

   } // if ( MPI_Rank == 0 )



// gravity solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )

// errors
// ------------------------------
   if ( OPT__EXT_ACC )
      Aux_Error( ERROR_INFO, "ELBDM does NOT support OPT__EXT_ACC --> use OPT__EXT_POT instead !!\n" );


// warnings
// ------------------------------
#  if ( GRA_GHOST_SIZE != 0  &&  !defined STORE_POT_GHOST )
#     warning : WARNING : GRA_GHOST_SIZE != 0 in ELBDM when STORE_POT_GHOST is off !!
#  endif


#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#endif // GRAVITY



// particle
// =======================================================================================
#ifdef PARTICLE

// errors
// ------------------------------
#  ifndef GRAVITY
#     error : ERROR : currently PARTICLE must work with GRAVITY !!
#  endif

#  ifdef COMOVING
#     error : ERROR : currently PARTICLE dost NOT support COMOVING !!
#  endif

#  if ( !defined SERIAL  &&  !defined LOAD_BALANCE )
#     error : ERROR : PARTICLE must work with either SERIAL or LOAD_BALANCE !!
#  endif

   if ( OPT__INIT != INIT_BY_RESTART )
   {
      if ( amr->Par->Init == PAR_INIT_BY_RESTART )    Aux_Error( ERROR_INFO, "PAR_INIT == RESTART but OPT__INIT != RESTART !!\n" );

      if ( amr->Par->NPar_Active_AllRank < 0 )
         Aux_Error( ERROR_INFO, "total number of particles in all MPI ranks = %ld < 0 !!\n",
                    amr->Par->NPar_Active_AllRank );

      if ( amr->Par->NPar_AcPlusInac < 0  ||  amr->Par->NPar_AcPlusInac > amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "incorrect total number of particles (%ld) in MPI rank %d !!\n",
                    amr->Par->NPar_AcPlusInac, MPI_Rank );
   }

#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "PAR_IMPROVE_ACC must work with STORE_POT_GHOST !!\n" );
#  endif

   if ( amr->Par->ImproveAcc  &&  amr->Par->Interp == 1 )
      Aux_Error( ERROR_INFO, "PAR_IMPROVE_ACC does NOT work with PAR_INTERP == 1 (NGP) !!\n" );

#  ifndef STORE_PAR_ACC
   if ( DT__PARACC != 0.0 )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) is NOT supported when STORE_PAR_ACC is off !!\n", DT__PARACC );
#  endif

   for (int d=0; d<3; d++)
   {
//    we have assumed that OPT__BC_FLU[2*d] == OPT__BC_FLU[2*d+1] when adopting the periodic BC
      if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC  &&  NX0_TOT[d]/PS2 == 1 )
         Aux_Error( ERROR_INFO, "\"%s\" does NOT work for NX0_TOT[%d] = 2*PATCH_SIZE when periodic BC is adopted !!\n",
                    "Par_MassAssignment()", d );
   }


// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( DT__PARVEL > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PARVEL (%13.7e) is not within the normal range [0.0~~1.0] !!\n", DT__PARVEL );

   if ( DT__PARACC > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PARACC (%13.7e) is not within the normal range [0.0~1.0] !!\n", DT__PARACC );

   if ( OPT__OVERLAP_MPI )
      Aux_Message( stderr, "WARNING : PARTICLE does not support OPT__OVERLAP_MPI !!\n" );

#  ifdef STORE_POT_GHOST
   if ( !amr->Par->ImproveAcc )
      Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless when PAR_IMPROVE_ACC is disabled !!\n" );
#  endif

   if ( OPT__GRA_P5_GRADIENT )
      Aux_Message( stderr, "WARNING : currently \"%s\" is not applied to particle update !!\n", "OPT__GRA_P5_GRADIENT" );

   } // if ( MPI_Rank == 0 )


#else // #ifdef PARTICLE


// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  ifdef STORE_POT_GHOST
   Aux_Message( stderr, "WARNING : currently STORE_POT_GHOST is useless when PARTICLE is disabled !!\n" );
#  endif

   }


#endif // PARTICLE



// Grackle
// =======================================================================================
#ifdef SUPPORT_GRACKLE

// errors
// ------------------------------
   /*
   if ( CHE_GPU_NPGROUP % GPU_NSTREAM != 0 )
      Aux_Error( ERROR_INFO, "CHE_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 CHE_GPU_NPGROUP, GPU_NSTREAM );
                 */

#  if ( EOS != EOS_GAMMA )
#     error : ERROR : SUPPORT_GRACKLE must work with EOS_GAMMA !!
#  endif

#  ifdef OPENMP
   if ( CHE_GPU_NPGROUP < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "CHE_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", CHE_GPU_NPGROUP, OMP_NTHREAD );
#  endif

// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( OPT__OVERLAP_MPI )
      Aux_Message( stderr, "WARNING : currently SUPPORT_GRACKLE does not support \"%s\" !!\n", "OPT__OVERLAP_MPI" );

   if ( GRACKLE_PRIMORDIAL > 0 )
      Aux_Message( stderr, "WARNING : adiabatic index gamma is currently fixed to %13.7e for Grackle !!\n", GAMMA );

   } // if ( MPI_Rank == 0 )

#endif // SUPPORT_GRACKLE



// source terms
// =======================================================================================

// errors
// ------------------------------
#  if ( SRC_GHOST_SIZE != 0 )
#     error : ERROR : SRC_GHOST_SIZE must be zero for now !!
#  endif

#  if ( MODEL != HYDRO )
   if ( SrcTerms.Deleptonization )
      Aux_Error( ERROR_INFO, "SRC_DELEPTONIZATION is only supported in HYDRO !!\n" );
#  endif

   if ( SRC_GPU_NPGROUP % GPU_NSTREAM != 0 )
      Aux_Error( ERROR_INFO, "SRC_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 SRC_GPU_NPGROUP, GPU_NSTREAM );

#  ifdef OPENMP
   if ( SRC_GPU_NPGROUP < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "SRC_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", SRC_GPU_NPGROUP, OMP_NTHREAD );
#  endif

// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

   } // if ( MPI_Rank == 0 )



// star formation
// =======================================================================================
#ifdef STAR_FORMATION

// errors
// ------------------------------
#  ifndef PARTICLE
#     error : STAR_FORMATION must work with PARTICLE !!
#  endif

#  if ( defined STORE_PAR_ACC  &&  !defined STORE_POT_GHOST )
#     error : STAR_FORMATION + STORE_PAR_ACC must work with STORE_POT_GHOST !!
#  endif

// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( SF_CREATE_STAR_MIN_LEVEL > MAX_LEVEL )
      Aux_Message( stderr, "WARNING : SF_CREATE_STAR_MIN_LEVEL (%d) > MAX_LEVEL (%d) --> no star particles will form !!\n",
                   SF_CREATE_STAR_MIN_LEVEL, MAX_LEVEL );

   if ( SF_CREATE_STAR_SCHEME == SF_CREATE_STAR_SCHEME_AGORA  &&  !SF_CREATE_STAR_DET_RANDOM )
   {
      Aux_Message( stderr, "WARNING : SF_CREATE_STAR_SCHEME == 1 will break bitwise reproducibility due to the \n" );
      Aux_Message( stderr, "          random values used for the stochastic star formation !!\n" );
      Aux_Message( stderr, "          --> Enable \"SF_CREATE_STAR_DET_RANDOM\" if reproducibility is of great concern\n" );
   }

   } // if ( MPI_Rank == 0 )

#endif // ifdef STAR_FORMATION


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ... done\n" );

} // FUNCTION : Aux_Check_Parameter
