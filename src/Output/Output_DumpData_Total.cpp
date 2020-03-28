#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif



/*======================================================================================================
Procedure for outputting new variables:
1. Add the new variables into one of the header sections
2. Add the corresponding (i) fread (ii) CompareVar in "Init_ByRestart_v2/Load_Parameter_After_2000"
======================================================================================================*/




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Total (FormatVersion = 2220)
// Description :  Output all simulation data in the binary form, which can be used as a restart file
//
// Note        :  1. This output format is deprecated and is mainly used for debugging only
//                   --> Use HDF5 format instead (OPT__OUTPUT_TOTAL = 1)
//
// Parameter   :  FileName : Name of the output file
//
// Revision    :  2110 : 2016/10/03 --> output HUBBLE0, OPT__UNIT, UNIT_L/M/T/V/D/E, MOLECULAR_WEIGHT
//                2120 : 2017/02/14 --> output passive grid and particle variables
//                2130 : 2017/08/09 --> output dTime_AllLv
//                2131 : 2017/12/05 --> no longer define INTEL
//                2200 : 2018/07/15 --> replace PAR_NVAR and PAR_NPASSIVE by PAR_NATT_STORED and PAR_NATT_USER;
//                                      use the new infrastructure for adding user-defined grid fields and
//                                      particle attributes
//                                      --> imcompatible with version 2131 for the data with user-defined particle
//                                          attributes as the order of their indices may be different
//                2201 : 2018/12/12 --> always set EP_COEFF=NULL_REAL since this variable no longer exists
//                2202 : 2018/12/15 --> set WAF-related variables to arbitrary values since they no longer exist
//                2203 : 2018/12/27 --> replace GRA_BLOCK_SIZE_Z by GRA_BLOCK_SIZE
//                2210 : 2019/06/07 --> support MHD
//                2220 : 2020/08/25 --> output EOS
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Total( const char *FileName )
{

// output data in the HDF5 format
#  ifdef SUPPORT_HDF5
   if ( OPT__OUTPUT_TOTAL == OUTPUT_FORMAT_HDF5 )
   {
      Output_DumpData_Total_HDF5( FileName );
      return;
   }
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the target file already exists
   if ( Aux_CheckFileExist(FileName)  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// get the total number of patches that have no son
   int NDataPatch_Local[NLEVEL] = { 0 };
   int NDataPatch_Total[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      if ( amr->patch[0][lv][PID]->son == -1 )  NDataPatch_Local[lv] ++;
   }

   MPI_Reduce( NDataPatch_Local, NDataPatch_Total, NLEVEL, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );


// get the number of partices in each rank and the corresponding global particle index offset
#  ifdef PARTICLE
   long GParID_Offset = 0;                   // GParID = global particle index (==> unique for each particle)
   long NPar_EachRank[MPI_NRank];

   MPI_Allgather( &amr->Par->NPar_Active, 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

   for (int r=0; r<MPI_Rank; r++)   GParID_Offset += NPar_EachRank[r];
#  endif


   FILE *File = NULL;
   long ExpectFileSize;
#  ifdef PARTICLE
   long FileOffset_Particle;
#  endif

   if ( MPI_Rank == 0 )
   {
      File = fopen( FileName, "wb" );

//    check the size of different data types
//    =================================================================================================
      const int size_bool   = sizeof( bool   );
      const int size_int    = sizeof( int    );
      const int size_uint   = sizeof( uint   );
      const int size_long   = sizeof( long   );
      const int size_ulong  = sizeof( ulong  );
      const int size_real   = sizeof( real   );
      const int size_double = sizeof( double );

      if ( size_int != size_uint )
         Aux_Error( ERROR_INFO, "sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );

      if ( size_long != size_ulong )
         Aux_Error( ERROR_INFO, "sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );


//    set the size of different headers (in bytes)
//    =================================================================================================
      const long HeaderSize_Format      =  256;
      const long HeaderSize_Makefile    =  256;
      const long HeaderSize_Constant    =  512;
      const long HeaderSize_Parameter   = 1024;
      const long HeaderSize_SimuInfo    = 1024;

      const long HeaderOffset_Format    = 0;    // it must be zero
      const long HeaderOffset_Makefile  = HeaderOffset_Format    + HeaderSize_Format;
      const long HeaderOffset_Constant  = HeaderOffset_Makefile  + HeaderSize_Makefile;
      const long HeaderOffset_Parameter = HeaderOffset_Constant  + HeaderSize_Constant;
      const long HeaderOffset_SimuInfo  = HeaderOffset_Parameter + HeaderSize_Parameter;

      const long HeaderSize_Total       = HeaderOffset_SimuInfo  + HeaderSize_SimuInfo;


//    fill in the entire header in advance (just for clarification)
//    =================================================================================================
      char *OutputBuf = new char [HeaderSize_Total];

      for (int t=0; t<HeaderSize_Total; t++)   OutputBuf[t] = 'b';   // 'b' is chosen arbitrarily ...

      fwrite( OutputBuf, sizeof(char), HeaderSize_Total, File );

      delete [] OutputBuf;


//    calculate the total file size
//    =================================================================================================
      int PatchDataSize, NGridVar=NCOMP_TOTAL;

#     ifdef GRAVITY
      if ( OPT__OUTPUT_POT )                                NGridVar ++;
#     endif
#     ifdef PARTICLE
      if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )   NGridVar ++;
#     endif
#     ifdef MHD
      if ( OPT__OUTPUT_CC_MAG )                             NGridVar += 3;
#     endif

      PatchDataSize  = CUBE(PS1)*NGridVar*sizeof(real);
#     ifdef MHD
      PatchDataSize += PS1P1*SQR(PS1)*NCOMP_MAG*sizeof(real);
#     endif
      ExpectFileSize = HeaderSize_Total;

      for (int lv=0; lv<NLEVEL; lv++)
      {
         ExpectFileSize += (long)NPatchTotal[lv]*4*sizeof(int);      // 4 = corner(3) + son(1)
         ExpectFileSize += (long)NDataPatch_Total[lv]*PatchDataSize;
      }

#     ifdef PARTICLE
      for (int lv=0; lv<NLEVEL; lv++)
         ExpectFileSize += (long)NDataPatch_Total[lv]*2*sizeof(long);   // 2 = NPar + starting particle index (leaf patches only)

      FileOffset_Particle = ExpectFileSize;  // file offset at the beginning of particle data

      ExpectFileSize += (long)PAR_NATT_STORED*amr->Par->NPar_Active_AllRank*sizeof(real);
#     endif


//    a. output the information of data format
//    =================================================================================================
      const long FormatVersion = 2220;
      const long CheckCode     = 123456789;

      fseek( File, HeaderOffset_Format, SEEK_SET );

      fwrite( &FormatVersion,             sizeof(long),                    1,             File );
      fwrite( &HeaderSize_Format,         sizeof(long),                    1,             File );
      fwrite( &HeaderSize_Makefile,       sizeof(long),                    1,             File );
      fwrite( &HeaderSize_Constant,       sizeof(long),                    1,             File );
      fwrite( &HeaderSize_Parameter,      sizeof(long),                    1,             File );
      fwrite( &HeaderSize_SimuInfo,       sizeof(long),                    1,             File );
      fwrite( &CheckCode,                 sizeof(long),                    1,             File );
      fwrite( &size_bool,                 sizeof(int),                     1,             File );
      fwrite( &size_int,                  sizeof(int),                     1,             File );
      fwrite( &size_long,                 sizeof(int),                     1,             File );
      fwrite( &size_real,                 sizeof(int),                     1,             File );
      fwrite( &size_double,               sizeof(int),                     1,             File );


//    b. output the simulation options and parameters defined in the Makefile
//    =================================================================================================
      if ( ftell(File) > HeaderOffset_Makefile )
         Aux_Error( ERROR_INFO, "Current file position (%ld) > correct offset (%ld) !!\n",
                    ftell(File), HeaderOffset_Makefile );

      fseek( File, HeaderOffset_Makefile, SEEK_SET );

#     ifdef MODEL
      const int  model               = MODEL;
#     else
      const int  model               = NULL_INT;
#     endif

#     ifdef GRAVITY
      const bool gravity             = true;
#     else
      const bool gravity             = false;
#     endif

#     ifdef POT_SCHEME
      const int  pot_scheme          = POT_SCHEME;
#     else
      const int  pot_scheme          = NULL_INT;
#     endif

      const bool individual_timestep = ( OPT__DT_LEVEL != DT_LEVEL_SHARED );

#     ifdef COMOVING
      const bool comoving            = true;
#     else
      const bool comoving            = false;
#     endif

#     ifdef FLU_SCHEME
      const int  flu_scheme          = FLU_SCHEME;
#     else
      const int  flu_scheme          = NULL_INT;
#     endif

#     ifdef LR_SCHEME
      const int  lr_scheme           = LR_SCHEME;
#     else
      const int  lr_scheme           = NULL_INT;
#     endif

#     ifdef RSOLVER
      const int  rsolver             = RSOLVER;
#     else
      const int  rsolver             = NULL_INT;
#     endif

#     ifdef GPU
      const bool gpu                 = true;
#     else
      const bool gpu                 = false;
#     endif

//    gamer_optimization is not used anymore after output version >= 2101
      const bool gamer_optimization  = NULL_BOOL;

#     ifdef GAMER_DEBUG
      const bool gamer_debug         = true;
#     else
      const bool gamer_debug         = false;
#     endif

#     ifdef TIMING
      const bool timing              = true;
#     else
      const bool timing              = false;
#     endif

#     ifdef TIMING_SOLVER
      const bool timing_solver       = true;
#     else
      const bool timing_solver       = false;
#     endif

//    compilation option "INTEL" has been removed
      const bool intel               = false;

#     ifdef FLOAT8
      const bool float8              = true;
#     else
      const bool float8              = false;
#     endif

#     ifdef SERIAL
      const bool serial              = true;
#     else
      const bool serial              = false;
#     endif

#     ifdef LOAD_BALANCE
      const int load_balance         = LOAD_BALANCE;
#     else
      const int load_balance         = 0;
#     endif

#     ifdef OVERLAP_MPI
      const bool overlap_mpi         = true;
#     else
      const bool overlap_mpi         = false;
#     endif

#     ifdef OPENMP
      const bool openmp              = true;
#     else
      const bool openmp              = false;
#     endif

#     ifdef GPU
      const int  gpu_arch            = GPU_ARCH;
#     else
      const int  gpu_arch            = NULL_INT;
#     endif

#     ifdef STORE_POT_GHOST
      const bool store_pot_ghost     = true;
#     else
      const bool store_pot_ghost     = false;
#     endif

#     ifdef UNSPLIT_GRAVITY
      const bool unsplit_gravity     = true;
#     else
      const bool unsplit_gravity     = false;
#     endif

#     ifdef PARTICLE
      const bool particle            = true;
#     else
      const bool particle            = false;
#     endif

      const int ncomp_passive        = NCOMP_PASSIVE;

#     ifdef CONSERVE_MASS
      const bool conserve_mass       = true;
#     else
      const bool conserve_mass       = false;
#     endif

#     ifdef LAPLACIAN_4TH
      const bool laplacian_4th       = true;
#     else
      const bool laplacian_4th       = false;
#     endif

#     ifdef QUARTIC_SELF_INTERACTION
      const bool self_interaction    = true;
#     else
      const bool self_interaction    = false;
#     endif

#     ifdef LAOHU
      const bool laohu               = true;
#     else
      const bool laohu               = false;
#     endif

#     ifdef SUPPORT_HDF5
      const bool support_hdf5        = true;
#     else
      const bool support_hdf5        = false;
#     endif

#     ifdef MHD
      const bool mhd                 = true;
#     else
      const bool mhd                 = false;
#     endif

#     if ( MODEL == HYDRO )
      const int eos                  = EOS;
#     else
      const int eos                  = NULL_INT;
#     endif

      const int nlevel               = NLEVEL;
      const int max_patch            = MAX_PATCH;

      fwrite( &model,                     sizeof(int),                     1,             File );
      fwrite( &gravity,                   sizeof(bool),                    1,             File );
      fwrite( &pot_scheme,                sizeof(int),                     1,             File );
      fwrite( &individual_timestep,       sizeof(bool),                    1,             File );
      fwrite( &comoving,                  sizeof(bool),                    1,             File );
      fwrite( &flu_scheme,                sizeof(int),                     1,             File );
      fwrite( &lr_scheme,                 sizeof(int),                     1,             File );
      fwrite( &rsolver,                   sizeof(int),                     1,             File );
      fwrite( &gpu,                       sizeof(bool),                    1,             File );
      fwrite( &gamer_optimization,        sizeof(bool),                    1,             File );
      fwrite( &gamer_debug,               sizeof(bool),                    1,             File );
      fwrite( &timing,                    sizeof(bool),                    1,             File );
      fwrite( &timing_solver,             sizeof(bool),                    1,             File );
      fwrite( &intel,                     sizeof(bool),                    1,             File );
      fwrite( &float8,                    sizeof(bool),                    1,             File );
      fwrite( &serial,                    sizeof(bool),                    1,             File );
      fwrite( &load_balance,              sizeof(int),                     1,             File );
      fwrite( &overlap_mpi,               sizeof(bool),                    1,             File );
      fwrite( &openmp,                    sizeof(bool),                    1,             File );
      fwrite( &gpu_arch,                  sizeof(int),                     1,             File );
      fwrite( &nlevel,                    sizeof(int),                     1,             File );
      fwrite( &max_patch,                 sizeof(int),                     1,             File );
      fwrite( &store_pot_ghost,           sizeof(bool),                    1,             File );
      fwrite( &unsplit_gravity,           sizeof(bool),                    1,             File );
      fwrite( &particle,                  sizeof(bool),                    1,             File );
      fwrite( &ncomp_passive,             sizeof(int),                     1,             File );
      fwrite( &conserve_mass,             sizeof(bool),                    1,             File );
      fwrite( &laplacian_4th,             sizeof(bool),                    1,             File );
      fwrite( &self_interaction,          sizeof(bool),                    1,             File );
      fwrite( &laohu,                     sizeof(bool),                    1,             File );
      fwrite( &support_hdf5,              sizeof(bool),                    1,             File );
      fwrite( &mhd,                       sizeof(bool),                    1,             File );
      fwrite( &eos,                       sizeof(int),                     1,             File );


//    c. output the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
//    =================================================================================================
      if ( ftell(File) > HeaderOffset_Constant )
         Aux_Error( ERROR_INFO, "Current file position (%ld) > correct offset (%ld) !!\n",
                    ftell(File), HeaderOffset_Constant );

      fseek( File, HeaderOffset_Constant, SEEK_SET );

      const int    ncomp_fluid           = NCOMP_FLUID;
      const int    patch_size            = PATCH_SIZE;
#     if ( MODEL == HYDRO )
      const double min_pres              = MIN_PRES;
#     else
      const double min_pres              = NULL_REAL;
#     endif
      const int    flu_ghost_size        = FLU_GHOST_SIZE;

#     ifdef GRAVITY
      const int    pot_ghost_size        = POT_GHOST_SIZE;
      const int    gra_ghost_size        = GRA_GHOST_SIZE;
#     else
      const int    pot_ghost_size        = NULL_INT;
      const int    gra_ghost_size        = NULL_INT;
#     endif

      const bool   enforce_positive      = true;

#     ifdef CHAR_RECONSTRUCTION
      const bool   char_reconstruction   = true;
#     else
      const bool   char_reconstruction   = false;
#     endif

#     ifdef CHECK_INTERMEDIATE
      const int    check_intermediate    = CHECK_INTERMEDIATE;
#     else
      const int    check_intermediate    = NULL_INT;
#     endif

#     ifdef HLL_NO_REF_STATE
      const bool   hll_no_ref_state      = true;
#     else
      const bool   hll_no_ref_state      = false;
#     endif

#     ifdef HLL_INCLUDE_ALL_WAVES
      const bool   hll_include_all_waves = true;
#     else
      const bool   hll_include_all_waves = false;
#     endif

      const bool   waf_dissipate_uesless = NULL_BOOL;    // this variable no longer exists
      const double max_error             = MAX_ERROR;
      const int    flu_block_size_x      = FLU_BLOCK_SIZE_X;
      const int    flu_block_size_y      = FLU_BLOCK_SIZE_Y;

#     ifdef USE_PSOLVER_10TO14
      const bool   use_psolver_10to14    = true;
#     else
      const bool   use_psolver_10to14    = false;
#     endif

#     ifdef POT_BLOCK_SIZE_X
      const int   pot_block_size_x       = POT_BLOCK_SIZE_X;
#     else
      const int   pot_block_size_x       = NULL_INT;
#     endif

#     ifdef POT_BLOCK_SIZE_Z
      const int   pot_block_size_z       = POT_BLOCK_SIZE_Z;
#     else
      const int   pot_block_size_z       = NULL_INT;
#     endif

#     ifdef GRA_BLOCK_SIZE
      const int   gra_block_size         = GRA_BLOCK_SIZE;
#     else
      const int   gra_block_size         = NULL_INT;
#     endif

#     ifdef PARTICLE
      const int    par_natt_stored       = PAR_NATT_STORED;
      const int    par_natt_user         = PAR_NATT_USER;
#     else
      const int    par_natt_stored       = NULL_INT;
      const int    par_natt_user         = NULL_INT;
#     endif

      fwrite( &ncomp_fluid,               sizeof(int),                     1,             File );
      fwrite( &patch_size,                sizeof(int),                     1,             File );
      fwrite( &min_pres,                  sizeof(double),                  1,             File );
      fwrite( &flu_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &pot_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &gra_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &enforce_positive,          sizeof(bool),                    1,             File );
      fwrite( &char_reconstruction,       sizeof(bool),                    1,             File );
      fwrite( &check_intermediate,        sizeof(int),                     1,             File );
      fwrite( &hll_no_ref_state,          sizeof(bool),                    1,             File );
      fwrite( &hll_include_all_waves,     sizeof(bool),                    1,             File );
      fwrite( &waf_dissipate_uesless,     sizeof(bool),                    1,             File );
      fwrite( &max_error,                 sizeof(double),                  1,             File );
      fwrite( &flu_block_size_x,          sizeof(int),                     1,             File );
      fwrite( &flu_block_size_y,          sizeof(int),                     1,             File );
      fwrite( &use_psolver_10to14,        sizeof(bool),                    1,             File );
      fwrite( &pot_block_size_x,          sizeof(int),                     1,             File );
      fwrite( &pot_block_size_z,          sizeof(int),                     1,             File );
      fwrite( &gra_block_size,            sizeof(int),                     1,             File );
      fwrite( &par_natt_stored,           sizeof(int),                     1,             File );
      fwrite( &par_natt_user,             sizeof(int),                     1,             File );


//    d. output the simulation parameters recorded in the file "Input__Parameter"
//    *** we always typecast "enum" to "int" when outputting data to ensure the data size ***
//    =================================================================================================
      if ( ftell(File) > HeaderOffset_Parameter )
         Aux_Error( ERROR_INFO, "Current file position (%ld) > correct offset (%ld) !!\n",
                    ftell(File), HeaderOffset_Parameter );

      fseek( File, HeaderOffset_Parameter, SEEK_SET );

      const int    mpi_nrank                 = MPI_NRank;
      const int    mpi_nrank_x[3]            = { MPI_NRank_X[0], MPI_NRank_X[1], MPI_NRank_X[2] };
      const int    opt__output_total         = (int)OPT__OUTPUT_TOTAL;
      const int    opt__output_part          = (int)OPT__OUTPUT_PART;
      const int    opt__output_mode          = (int)OPT__OUTPUT_MODE;
      const int    opt__flu_int_scheme       = (int)OPT__FLU_INT_SCHEME;
      const int    opt__ref_flu_int_scheme   = (int)OPT__REF_FLU_INT_SCHEME;

#     ifndef COMOVING
      const double OMEGA_M0                  = NULL_REAL;
      const double HUBBLE0                   = NULL_REAL;
      const double DT__MAX_DELTA_A           = NULL_REAL;
#     endif

#     ifdef GRAVITY
      const int    opt__pot_int_scheme       = (int)OPT__POT_INT_SCHEME;
      const int    opt__rho_int_scheme       = (int)OPT__RHO_INT_SCHEME;
      const int    opt__gra_int_scheme       = (int)OPT__GRA_INT_SCHEME;
      const int    opt__ref_pot_int_scheme   = (int)OPT__REF_POT_INT_SCHEME;
#     else
      const double DT__GRAVITY               = NULL_REAL;
      const double NEWTON_G                  = NULL_REAL;
      const int    POT_GPU_NPGROUP           = NULL_INT;
      const bool   OPT__OUTPUT_POT           = false;
      const bool   OPT__GRA_P5_GRADIENT      = NULL_BOOL;
      const double SOR_OMEGA                 = NULL_REAL;
      const int    SOR_MAX_ITER              = NULL_INT;
      const int    SOR_MIN_ITER              = NULL_INT;
      const double MG_TOLERATED_ERROR        = NULL_REAL;
      const int    MG_MAX_ITER               = NULL_INT;
      const int    MG_NPRE_SMOOTH            = NULL_INT;
      const int    MG_NPOST_SMOOTH           = NULL_INT;
      const int    opt__pot_int_scheme       = NULL_INT;
      const int    opt__rho_int_scheme       = NULL_INT;
      const int    opt__gra_int_scheme       = NULL_INT;
      const int    opt__ref_pot_int_scheme   = NULL_INT;
#     endif // #ifdef GRAVITY

#     ifdef LOAD_BALANCE
      const double lb_wli_max                = amr->LB->WLI_Max;
#     else
      const double lb_wli_max                = NULL_REAL;
#     endif

      const double EP_COEFF_useless          = NULL_REAL;   // this variable no longer exists
      const int    opt__waf_limiter_useless  = NULL_INT;    // this variable no longer exists

#     if ( MODEL == HYDRO )
      const int    opt__lr_limiter           = (int)OPT__LR_LIMITER;

//    convert OPT__1ST_FLUX_CORR to bool to be consistent with the old format where OPT__1ST_FLUX_CORR is bool instead of int
      const bool   opt__1st_flux_corr        = (bool)OPT__1ST_FLUX_CORR;
      const int    opt__1st_flux_corr_scheme = (int)OPT__1ST_FLUX_CORR_SCHEME;
#     else
      const bool   OPT__FLAG_PRES_GRADIENT   = NULL_BOOL;
      const double GAMMA                     = NULL_REAL;
      const double MOLECULAR_WEIGHT          = NULL_REAL;
      const double MINMOD_COEFF              = NULL_REAL;
      const int    opt__lr_limiter           = NULL_INT;
      const bool   opt__1st_flux_corr        = NULL_BOOL;
      const int    opt__1st_flux_corr_scheme = NULL_INT;
#     endif

#     if ( MODEL != ELBDM )
      const double DT__PHASE                 = NULL_REAL;
      const bool   OPT__FLAG_ENGY_DENSITY    = NULL_BOOL;
      const bool   OPT__INT_PHASE            = NULL_BOOL;
      const double ELBDM_MASS                = NULL_REAL;
      const double ELBDM_PLANCK_CONST        = NULL_REAL;
#     endif

#     ifdef PARTICLE
      const int    opt__output_par_dens      = (int)OPT__OUTPUT_PAR_DENS;
#     else
      const int    opt__output_par_dens      = 0;
#     endif

#     ifndef MHD
      const bool   OPT__OUTPUT_CC_MAG        = false;
      const double UNIT_B                    = NULL_REAL;
#     endif

      const bool   dummy_bool                = false;

      fwrite( &BOX_SIZE,                  sizeof(double),                  1,             File );
      fwrite(  NX0_TOT,                   sizeof(int),                     3,             File );
      fwrite( &mpi_nrank,                 sizeof(int),                     1,             File );
      fwrite(  mpi_nrank_x,               sizeof(int),                     3,             File );
      fwrite( &OMP_NTHREAD,               sizeof(int),                     1,             File );
      fwrite( &END_T,                     sizeof(double),                  1,             File );
      fwrite( &END_STEP,                  sizeof(long),                    1,             File );
      fwrite( &OMEGA_M0,                  sizeof(double),                  1,             File );
      fwrite( &DT__FLUID,                 sizeof(double),                  1,             File );
      fwrite( &DT__GRAVITY,               sizeof(double),                  1,             File );
      fwrite( &DT__PHASE,                 sizeof(double),                  1,             File );
      fwrite( &DT__MAX_DELTA_A,           sizeof(double),                  1,             File );
      fwrite( &dummy_bool,                sizeof(bool),                    1,             File );
      fwrite( &OPT__DT_USER,              sizeof(bool),                    1,             File );
      fwrite( &REGRID_COUNT,              sizeof(int),                     1,             File );
      fwrite( &FLAG_BUFFER_SIZE,          sizeof(int),                     1,             File );
      fwrite( &MAX_LEVEL,                 sizeof(int),                     1,             File );
      fwrite( &OPT__FLAG_RHO,             sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_RHO_GRADIENT,    sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_PRES_GRADIENT,   sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_ENGY_DENSITY,    sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_USER,            sizeof(bool),                    1,             File );
      fwrite( &lb_wli_max,                sizeof(double),                  1,             File );
      fwrite( &GAMMA,                     sizeof(double),                  1,             File );
      fwrite( &MINMOD_COEFF,              sizeof(double),                  1,             File );
      fwrite( &EP_COEFF_useless,          sizeof(double),                  1,             File );
      fwrite( &opt__lr_limiter,           sizeof(int),                     1,             File );
      fwrite( &opt__waf_limiter_useless,  sizeof(int),                     1,             File );
      fwrite( &ELBDM_MASS,                sizeof(double),                  1,             File );
      fwrite( &ELBDM_PLANCK_CONST,        sizeof(double),                  1,             File );
      fwrite( &FLU_GPU_NPGROUP,           sizeof(int),                     1,             File );
      fwrite( &GPU_NSTREAM,               sizeof(int),                     1,             File );
      fwrite( &OPT__FIXUP_FLUX,           sizeof(bool),                    1,             File );
      fwrite( &OPT__FIXUP_RESTRICT,       sizeof(bool),                    1,             File );
      fwrite( &OPT__OVERLAP_MPI,          sizeof(bool),                    1,             File );
      fwrite( &NEWTON_G,                  sizeof(double),                  1,             File );
      fwrite( &SOR_OMEGA,                 sizeof(double),                  1,             File );
      fwrite( &SOR_MAX_ITER,              sizeof(int),                     1,             File );
      fwrite( &SOR_MIN_ITER,              sizeof(int),                     1,             File );
      fwrite( &MG_MAX_ITER,               sizeof(int),                     1,             File );
      fwrite( &MG_NPRE_SMOOTH,            sizeof(int),                     1,             File );
      fwrite( &MG_NPOST_SMOOTH,           sizeof(int),                     1,             File );
      fwrite( &MG_TOLERATED_ERROR,        sizeof(double),                  1,             File );
      fwrite( &POT_GPU_NPGROUP,           sizeof(int),                     1,             File );
      fwrite( &OPT__GRA_P5_GRADIENT,      sizeof(bool),                    1,             File );
      fwrite( &OPT__INT_TIME,             sizeof(bool),                    1,             File );
      fwrite( &OPT__INT_PHASE,            sizeof(bool),                    1,             File );
      fwrite( &opt__flu_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__pot_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__rho_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__gra_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__ref_flu_int_scheme,   sizeof(int),                     1,             File );
      fwrite( &opt__ref_pot_int_scheme,   sizeof(int),                     1,             File );
      fwrite( &opt__output_total,         sizeof(int),                     1,             File );
      fwrite( &opt__output_part,          sizeof(int),                     1,             File );
      fwrite( &OPT__OUTPUT_USER,          sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_BASE,          sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_POT,           sizeof(bool),                    1,             File );
      fwrite( &opt__output_mode,          sizeof(int),                     1,             File );
      fwrite( &OUTPUT_STEP,               sizeof(int),                     1,             File );
      fwrite( &OUTPUT_DT,                 sizeof(double),                  1,             File );
      fwrite( &OUTPUT_PART_X,             sizeof(double),                  1,             File );
      fwrite( &OUTPUT_PART_Y,             sizeof(double),                  1,             File );
      fwrite( &OUTPUT_PART_Z,             sizeof(double),                  1,             File );
      fwrite( &OPT__TIMING_BALANCE,       sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_BASEPS,        sizeof(bool),                    1,             File );
      fwrite( &opt__1st_flux_corr,        sizeof(bool),                    1,             File );
      fwrite( &opt__1st_flux_corr_scheme, sizeof(int),                     1,             File );
      fwrite( &opt__output_par_dens,      sizeof(int),                     1,             File );
      fwrite( &HUBBLE0,                   sizeof(double),                  1,             File );
      fwrite( &OPT__UNIT,                 sizeof(bool),                    1,             File );
      fwrite( &UNIT_L,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_M,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_T,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_V,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_D,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_E,                    sizeof(double),                  1,             File );
      fwrite( &UNIT_P,                    sizeof(double),                  1,             File );
      fwrite( &MOLECULAR_WEIGHT,          sizeof(double),                  1,             File );
      fwrite( &OPT__OUTPUT_CC_MAG,        sizeof(bool),                    1,             File );
      fwrite( &UNIT_B,                    sizeof(double),                  1,             File );


//    e. output the simulation information
//    =================================================================================================
      if ( ftell(File) > HeaderOffset_SimuInfo )
         Aux_Error( ERROR_INFO, "Current file position (%ld) > correct offset (%ld) !!\n",
                    ftell(File), HeaderOffset_SimuInfo );

      fseek( File, HeaderOffset_SimuInfo, SEEK_SET );

#     ifndef GRAVITY
      const double AveDensity_Init     = NULL_REAL;
#     endif
#     ifdef PARTICLE
      const long   NParAllRank         = amr->Par->NPar_Active_AllRank;
#     else
      const long   NParAllRank         = NULL_INT;
      const long   FileOffset_Particle = NULL_INT;
#     endif

      fwrite( &CheckCode,                 sizeof(long),                    1,             File );
      fwrite( &DumpID,                    sizeof(int),                     1,             File );
      fwrite( Time,                       sizeof(double),             NLEVEL,             File );
      fwrite( &Step,                      sizeof(long),                    1,             File );
      fwrite( NPatchTotal,                sizeof(int),                NLEVEL,             File );
      fwrite( NDataPatch_Total,           sizeof(int),                NLEVEL,             File );
      fwrite( AdvanceCounter,             sizeof(long),               NLEVEL,             File );
      fwrite( &AveDensity_Init,           sizeof(double),                  1,             File );
      fwrite( &NParAllRank,               sizeof(long),                    1,             File );
      fwrite( &FileOffset_Particle,       sizeof(long),                    1,             File );
      fwrite( dTime_AllLv,                sizeof(double),             NLEVEL,             File );


//    move the file position indicator to the end of the header ==> prepare to output patch data
      if ( ftell(File) > HeaderSize_Total )
         Aux_Error( ERROR_INFO, "Current file position (%ld) > correct offset (%ld) !!\n",
                    ftell(File), HeaderSize_Total );

      fseek( File, HeaderSize_Total, SEEK_SET );

      fclose( File );

   } // if ( MPI_Rank == 0 )


// f. output the simulation grid data
// =================================================================================================
   int Cr_and_Son[4];

#  ifdef PARTICLE
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;
   const bool TimingSendPar_No  = false;
   const bool PredictParPos_No  = false;   // particles synchronization is done in "Flu_CorrAfterAllSync()"
   const bool JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool SibBufPatch       = true;
   const bool FaSibBufPatch     = true;
#  else
   const bool SibBufPatch       = NULL_BOOL;
   const bool FaSibBufPatch     = NULL_BOOL;
#  endif
   real (*ParDensArray)[ CUBE(PS1) ] = ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_NONE ) ? NULL : new real [8][ CUBE(PS1) ];

   long NPar_and_GParID[2], GParID=GParID_Offset;
#  endif // #ifdef PARTICLE

#  ifdef MHD
   real (*CCMag)[ CUBE(PS1) ] = ( OPT__OUTPUT_CC_MAG ) ? new real [NCOMP_MAG][ CUBE(PS1) ] : NULL;
#  endif


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    f-1. initialize the particle density array (rho_ext) and collect particles from higher levels for outputting particle density
#     ifdef MASSIVE_PARTICLES
      if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
      {
         Prepare_PatchData_InitParticleDensityArray( lv );

         Par_CollectParticle2OneLevel( lv, PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                       TimingSendPar_No );
      }
#     endif

      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "ab" );

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
//             f1. prepare the particle density data on grids (only if there are leaf patches in this patch group)
#              ifdef MASSIVE_PARTICLES 
               if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE  &&  PID%8 == 0 )
               {
                  for (int PID_CheckSon=PID; PID_CheckSon<PID+8; PID_CheckSon++)
                  {
                     if ( amr->patch[0][lv][PID_CheckSon]->son == -1 )
                     {
//                      we do not check minimum density here (just because it's unnecessary)
                        Prepare_PatchData( lv, Time[lv], ParDensArray[0], NULL, 0, 1, &PID,
                                           ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_PAR_ONLY ) ? _PAR_DENS : _TOTAL_DENS,
                                           _NONE, OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No,
                                           OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );
                        break;
                     }
                  }
               }
#              endif // #ifdef MASSIVE_PARTICLES 


//             f2. prepare the cell-centered magnetic field for leaf patches
#              ifdef MHD
               if ( OPT__OUTPUT_CC_MAG  &&  amr->patch[0][lv][PID]->son == -1 )
               {
                  real CCMag_1Cell[NCOMP_MAG];
                  int  idx=0;

                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     MHD_GetCellCenteredBFieldInPatch( CCMag_1Cell, lv, PID, i, j, k, amr->MagSg[lv] );

                     for (int v=0; v<NCOMP_MAG; v++)  CCMag[v][idx] = CCMag_1Cell[v];

                     idx ++;
                  }
               }
#              endif


//             f3. output patch information
//             (father <-> son information will be re-constructed during the restart)
               Cr_and_Son[0] = amr->patch[0][lv][PID]->corner[0];
               Cr_and_Son[1] = amr->patch[0][lv][PID]->corner[1];
               Cr_and_Son[2] = amr->patch[0][lv][PID]->corner[2];
               Cr_and_Son[3] = amr->patch[0][lv][PID]->son;

               fwrite( Cr_and_Son, sizeof(int), 4, File );

#              ifdef PARTICLE
//             output particle information for leaf patches only
               if ( amr->patch[0][lv][PID]->son == -1 )
               {
                  NPar_and_GParID[0] = amr->patch[0][lv][PID]->NPar;
                  NPar_and_GParID[1] = GParID;

                  fwrite( NPar_and_GParID, sizeof(long), 2, File );

                  GParID += amr->patch[0][lv][PID]->NPar;
               }
#              endif


//             f4. output patch data only if it has no son
               if ( amr->patch[0][lv][PID]->son == -1 )
               {
//                f4-1. output fluid variables
                  fwrite( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid,    sizeof(real), CUBE(PS1)*NCOMP_TOTAL,    File );

//                f4-2. output gravitational potential
#                 ifdef GRAVITY
                  if ( OPT__OUTPUT_POT )
                  fwrite( amr->patch[ amr->PotSg[lv] ][lv][PID]->pot,      sizeof(real), CUBE(PS1),                File );
#                 endif

//                f4-3. output particle density depostied onto grids
#                 ifdef PARTICLE
                  if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
                  fwrite( ParDensArray[ PID%8 ],                           sizeof(real), CUBE(PS1),                File );
#                 endif

//                f4-4. output magnetic field
#                 ifdef MHD
//                cell-centered
                  if ( OPT__OUTPUT_CC_MAG )
                  fwrite( CCMag,                                           sizeof(real), CUBE(PS1)*NCOMP_MAG,      File );

//                face-centered
                  fwrite( amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic, sizeof(real), PS1P1*SQR(PS1)*NCOMP_MAG, File );
#                 endif
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

            fclose( File );


//          f5. free memory used for outputting particle density
#           ifdef MASSIVE_PARTICLES
            if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
            {
               Prepare_PatchData_FreeParticleDensityArray( lv );

               Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );
            }
#           endif
         } // if ( MPI_Rank == TargetMPIRank )

         MPI_Barrier( MPI_COMM_WORLD );

      } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)

// free memory
#  ifdef PARTICLE
   delete [] ParDensArray;
#  endif

#  ifdef MHD
   delete [] CCMag;
#  endif


// g. output particles
// =================================================================================================
#  ifdef PARTICLE
// check if all active particles are captured
   if ( GParID != GParID_Offset + amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "final GParID (%ld) != expect (%ld) !!\n",
                 GParID, GParID_Offset + amr->Par->NPar_Active );


// allocate I/O buffer (just for better I/O performance)
   const long ParBufSize = 10000000;   // number of particles dumped at a time

   real *ParBuf = new real [ParBufSize];


// set the file offset of particle data for each rank
   MPI_Bcast( &FileOffset_Particle, 1, MPI_LONG, 0, MPI_COMM_WORLD );


// output particle data (one attribute at a time to avoid creating holes in the file)
   const long ParDataSize1v = amr->Par->NPar_Active_AllRank*sizeof(real);

   long NParInBuf, ParID, FileOffset_ThisVar;
   int  NParThisPatch;

   for (int v=0; v<PAR_NATT_STORED; v++)
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         File = fopen( FileName, "ab" );

//       set file position indicator to the end of the current file and check whether it's consistent with expectation
         FileOffset_ThisVar = FileOffset_Particle + ParDataSize1v*v + GParID_Offset*sizeof(real);
         NParInBuf          = 0;

         fseek( File, 0, SEEK_END );

         if ( ftell(File) != FileOffset_ThisVar )
            Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                       FileName, ftell(File), FileOffset_ThisVar );

         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            NParThisPatch = amr->patch[0][lv][PID]->NPar;

//          check if the particle I/O buffer is exceeded (possible only if the number of particles in this patch > ParBufSize)
            if ( NParInBuf + NParThisPatch > ParBufSize )
            {
               Aux_Message( stderr, "ERROR : NParInBuf (%ld) + NParThisPatch (%d) = %ld > ParBufSize (%ld) ...\n",
                            NParInBuf, NParThisPatch, NParInBuf+NParThisPatch, ParBufSize );
               Aux_Message( stderr, "        ==> Please increase ParBufSize in the file \"%s\" !!\n", __FILE__ );
               MPI_Exit();
            }

//          store particle data into the I/O buffer
            for (int p=0; p<NParThisPatch; p++)
            {
               ParID = amr->patch[0][lv][PID]->ParList[p];

               ParBuf[ NParInBuf ++ ] = amr->Par->Attribute[v][ParID];
            }

//          store particle data from I/O buffer to disk
            if ( PID+1 == amr->NPatchComma[lv][1]  ||  NParInBuf + amr->patch[0][lv][PID+1]->NPar > ParBufSize )
            {
               fwrite( ParBuf, sizeof(real), NParInBuf, File );

               NParInBuf = 0;
            }
         } // for PID, lv

         fclose( File );
      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++), for (int v=0; v<PAR_NATT_STORED; v++)

   delete [] ParBuf;
#  endif // #ifdef PARTICLE


// check the file size
   if ( MPI_Rank == 0 )
   {
      FILE *FileCheck = fopen( FileName, "rb" );
      long  FileSize;

      fseek( FileCheck, 0, SEEK_END );
      FileSize = ftell( FileCheck );

      if ( FileSize != ExpectFileSize )
         Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                    FileName, FileSize, ExpectFileSize );

      fclose( FileCheck );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Total


