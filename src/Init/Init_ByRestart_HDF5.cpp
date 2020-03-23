#ifdef SUPPORT_HDF5

#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <typeinfo>

void FillIn_Makefile (  Makefile_t &Makefile  );
void FillIn_SymConst (  SymConst_t &SymConst  );
void FillIn_InputPara( InputPara_t &InputPara );

template <typename T>
static herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                         const T *ComprPtr, const int NCompr, const bool Fatal_Compr );
static void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive,
                          const int *SonList, const int (*CrList)[3],
                          const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                          const hid_t *H5_SetID_FCMag, const hid_t *H5_SpaceID_FCMag, const hid_t *H5_MemID_FCMag,
                          const int *NParList, real **ParBuf, long *NewParList, const hid_t *H5_SetID_ParData,
                          const hid_t H5_SpaceID_ParData, const long *GParID_Offset, const long NParThisRank );
static void Check_Makefile ( const char *FileName, const int FormatVersion );
static void Check_SymConst ( const char *FileName, const int FormatVersion );
static void Check_InputPara( const char *FileName, const int FormatVersion );
static void ResetParameter( const char *FileName, double *EndT, long *EndStep );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByRestart_HDF5
// Description :  Reload a previous HDF5 output as the initial condition
//
// Note        :  1. This function will be invoked by "Init_ByRestart" automatically if the restart file
//                   is in the HDF5 format
//                2. Only work for format version >= 2100 (PARTICLE only works for version >= 2200)
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
void Init_ByRestart_HDF5( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( MPI_Rank == 0 )
   {
      if ( !Aux_CheckFileExist(FileName) )
         Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" does not exist !!\n", FileName );

      if ( !H5Fis_hdf5(FileName) )
         Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" is not in the HDF5 format !!\n", FileName );
   }

#  if ( defined PARTICLE  &&  !defined SERIAL  &&  !defined LOAD_BALANCE )
#     error : ERROR : PARTICLE must work with either SERIAL or LOAD_BALANCE !!
#  endif

   MPI_Barrier( MPI_COMM_WORLD );


// 1. load the simulation info
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal             = true;
   const bool NonFatal             = false;
   const int  Model                = MODEL;
   const int  NCompFluid           = NCOMP_FLUID;
   const int  NCompPassive         = NCOMP_PASSIVE;
   const int  PatchSize            = PS1;
#  ifdef GRAVITY
   const int  Gravity              = 1;
#  else
   const int  Gravity              = 0;
#  endif
#  ifdef PARTICLE
   const int  Particle             = 1;
   const int  Par_NAttStored       = PAR_NATT_STORED;
#  else
   const int  Particle             = 0;
#  endif
#  if ( MODEL == HYDRO )
#  ifdef MHD
   const int  Magnetohydrodynamics = 1;
#  else
   const int  Magnetohydrodynamics = 0;
#  endif
#  ifdef COSMIC_RAY
   const int  CosmicRay            = 1;
#  else
   const int  CosmicRay            = 0;
#  endif
#  endif // #if ( MODEL == HYDRO )

   KeyInfo_t KeyInfo;

   hid_t  H5_FileID, H5_SetID_KeyInfo, H5_TypeID_KeyInfo, H5_SetID_Cr;
#  ifdef LOAD_BALANCE
   hid_t  H5_SetID_LBIdx;
#  else
   hid_t  H5_SetID_Son;
#  endif
   herr_t H5_Status;
   int    NLvRescale, NPatchAllLv, GID_LvStart[NLEVEL];
   int   *NullPtr = NULL;


// 1-1. open the HDF5 file
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 1-2. load the KeyInfo dataset and datatype
   H5_SetID_KeyInfo  = H5Dopen( H5_FileID, "Info/KeyInfo", H5P_DEFAULT );
   if ( H5_SetID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/KeyInfo" );

   H5_TypeID_KeyInfo = H5Dget_type( H5_SetID_KeyInfo );
   if ( H5_TypeID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "Info/KeyInfo" );


// 1-3. load all target fields in KeyInfo one-by-one (by all ranks)
   LoadField( "FormatVersion",  &KeyInfo.FormatVersion,  H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );

// format version for HDF5 output
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "      The format version of the HDF5 RESTART file = %ld\n", KeyInfo.FormatVersion );

      if ( KeyInfo.FormatVersion < 2100 )
         Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2100) !!\n" );

#     ifdef PARTICLE
      if ( KeyInfo.FormatVersion < 2200 )
         Aux_Error( ERROR_INFO, "unsupported data format version for PARTICLE (only support version >= 2200) !!\n" );
#     endif

      if ( KeyInfo.FormatVersion < 2300 )
         Aux_Message( stderr, "WARNING : loading user-defined fields or particle attributes from version < 2300 will likely fail !!\n" );

#     ifdef MHD
      if ( KeyInfo.FormatVersion < 2400 )
         Aux_Error( ERROR_INFO, "unsupported data format version for MHD (only support version >= 2400) !!\n" );
#     endif
   }

   MPI_Barrier( MPI_COMM_WORLD );

   LoadField( "Model",          &KeyInfo.Model,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &Model,         1,    Fatal );
   LoadField( "Gravity",        &KeyInfo.Gravity,        H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &Gravity,       1,    Fatal );
   LoadField( "Particle",       &KeyInfo.Particle,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &Particle,      1,    Fatal );
   LoadField( "NLevel",         &KeyInfo.NLevel,         H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,      -1, NonFatal );
   LoadField( "NCompFluid",     &KeyInfo.NCompFluid,     H5_SetID_KeyInfo, H5_TypeID_KeyInfo, NonFatal, &NCompFluid,    1,    Fatal );
   LoadField( "NCompPassive",   &KeyInfo.NCompPassive,   H5_SetID_KeyInfo, H5_TypeID_KeyInfo, NonFatal, &NCompPassive,  1,    Fatal );
   LoadField( "PatchSize",      &KeyInfo.PatchSize,      H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &PatchSize,     1,    Fatal );

// runtime NLEVEL must be >= loaded NLEVEL
   if      ( KeyInfo.NLevel > NLEVEL )
      Aux_Error( ERROR_INFO, "%s : RESTART file (%d) > runtime (%d) !!\n",
                 "NLEVEL", KeyInfo.NLevel, NLEVEL );
   else
   {
      NLvRescale = 1 << ( NLEVEL - KeyInfo.NLevel );

      if ( MPI_Rank == 0  &&  NLvRescale != 1 )
         Aux_Message( stderr, "WARNING : the NLEVEL rescale factor is set to %d\n", NLvRescale );
   }

   MPI_Barrier( MPI_COMM_WORLD );

   LoadField( "DumpID",               &KeyInfo.DumpID,               H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "NX0",                   KeyInfo.NX0,                  H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NX0_TOT,               3,    Fatal );
   LoadField( "BoxScale",              KeyInfo.BoxScale,             H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "NPatch",                KeyInfo.NPatch,               H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "CellScale",             KeyInfo.CellScale,            H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
#  if ( MODEL == HYDRO )
   if ( KeyInfo.FormatVersion >= 2400 )
   LoadField( "Magnetohydrodynamics", &KeyInfo.Magnetohydrodynamics, H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &Magnetohydrodynamics,  1,    Fatal );

   if ( KeyInfo.FormatVersion >= 2421 )
   LoadField( "CosmicRay",            &KeyInfo.CosmicRay,            H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal, &CosmicRay,             1,    Fatal );
#  endif

   LoadField( "Step",                 &KeyInfo.Step,                 H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "AdvanceCounter",        KeyInfo.AdvanceCounter,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Par_NPar",             &KeyInfo.Par_NPar,             H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   if ( KeyInfo.FormatVersion >= 2300 )
   LoadField( "Par_NAttStored",       &KeyInfo.Par_NAttStored,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo, NonFatal, &Par_NAttStored,        1,    Fatal );
   else
   LoadField( "Par_NAttStored",       &KeyInfo.Par_NAttStored,       H5_SetID_KeyInfo, H5_TypeID_KeyInfo, NonFatal, &Par_NAttStored,        1, NonFatal );
#  endif

   LoadField( "BoxSize",               KeyInfo.BoxSize,              H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  amr->BoxSize,          3,    Fatal );
   LoadField( "Time",                  KeyInfo.Time,                 H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "CellSize",              KeyInfo.CellSize,             H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "dTime_AllLv",           KeyInfo.dTime_AllLv,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo, NonFatal,  NullPtr,              -1, NonFatal );
#  ifdef GRAVITY
   LoadField( "AveDens_Init",         &KeyInfo.AveDens_Init,         H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
#  endif

   LoadField( "CodeVersion",          &KeyInfo.CodeVersion,          H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );
   LoadField( "DumpWallTime",         &KeyInfo.DumpWallTime,         H5_SetID_KeyInfo, H5_TypeID_KeyInfo,    Fatal,  NullPtr,              -1, NonFatal );


// 1-4. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Fclose( H5_FileID );


// 1-5. set internal parameters
// 1-5-1. parameters must be reset
   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
      NPatchTotal[lv] = KeyInfo.NPatch[lv];
   }

#  ifdef PARTICLE
   amr->Par->NPar_Active_AllRank = KeyInfo.Par_NPar;
#  endif

// 1-5-2. parameters reset only when OPT__RESTART_RESET is disabled
   if ( ! OPT__RESTART_RESET )
   {
      for (int lv=0; lv<KeyInfo.NLevel; lv++)
      {
         Time          [lv] = KeyInfo.Time          [lv];
         AdvanceCounter[lv] = KeyInfo.AdvanceCounter[lv];

         if ( KeyInfo.FormatVersion >= 2250 )
         dTime_AllLv   [lv] = KeyInfo.dTime_AllLv   [lv];
      }

      Step            = KeyInfo.Step;
#     ifdef GRAVITY
      AveDensity_Init = KeyInfo.AveDens_Init;
#     endif
   }


// 1-6. set parameters in levels that do not exist in the input file
// --> assuming dTime_AllLv[] has been initialized as 0.0 properly
   for (int lv=KeyInfo.NLevel; lv<NLEVEL; lv++)
   {
      Time          [lv] = 0.0;
      NPatchTotal   [lv] = 0;
      AdvanceCounter[lv] = 0;
   }


// 1-7. set SgTime
   for (int lv=0; lv<NLEVEL; lv++)
   {
      amr->FluSgTime[lv][ amr->FluSg[lv] ] = Time[lv];
#     ifdef MHD
      amr->MagSgTime[lv][ amr->MagSg[lv] ] = Time[lv];
#     endif
#     ifdef GRAVITY
      amr->PotSgTime[lv][ amr->PotSg[lv] ] = Time[lv];
#     endif
   }


// 1-8. set the next dump ID
   if ( INIT_DUMPID < 0 )
      DumpID = ( OPT__RESTART_RESET ) ? 0 : KeyInfo.DumpID + 1;
   else
      DumpID = INIT_DUMPID;


// 1-9. reset parameters from the restart file
   if ( ! OPT__RESTART_RESET )   ResetParameter( FileName, &END_T, &END_STEP );


// 1-10. check all other simulation information (by rank 0 only)
   if ( MPI_Rank == 0 )
   {
      Check_Makefile ( FileName, KeyInfo.FormatVersion );
      Check_SymConst ( FileName, KeyInfo.FormatVersion );
      Check_InputPara( FileName, KeyInfo.FormatVersion );
   }


// 1-11. set the GID offset at different levels
   NPatchAllLv = 0;
   for (int lv=0; lv<NLEVEL; lv++)
   {
      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }

   MPI_Barrier( MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 2. load the tree information (load-balance indices, corner, son, ... etc) of all patches (by all ranks)
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

// 2-1. corner
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading corner table ...\n" );

// allocate memory
   int (*CrList_AllLv)[3] = new int [ NPatchAllLv ][3];

// load data
   H5_SetID_Cr = H5Dopen( H5_FileID, "Tree/Corner", H5P_DEFAULT );
   if ( H5_SetID_Cr < 0 )     Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Corner" );
   H5_Status = H5Dread( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Cr );

// rescale the loaded corner (necessary when KeyInfo.NLevel != NLEVEL)
   if ( NLvRescale != 1 )
   for (int GID=0; GID<NPatchAllLv; GID++)
   for (int d=0; d<3; d++)
      CrList_AllLv[GID][d] *= NLvRescale;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading corner table ... done\n" );


// 2-2. LBIdx
#  ifdef LOAD_BALANCE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading load-balance index table ...\n" );

   long *LBIdxList_AllLv;
   long *LBIdxList_EachLv         [NLEVEL];
   int  *LBIdxList_EachLv_IdxTable[NLEVEL];
   int   LoadIdx_Start[NLEVEL], LoadIdx_Stop[NLEVEL];

// 2-2-1. allocate memory
   LBIdxList_AllLv = new long [ NPatchAllLv ];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( lv < KeyInfo.NLevel )
      {
         LBIdxList_EachLv         [lv] = LBIdxList_AllLv + GID_LvStart[lv];
         LBIdxList_EachLv_IdxTable[lv] = new int [ NPatchTotal[lv] ];
      }

      else
      {
         LBIdxList_EachLv         [lv] = NULL;
         LBIdxList_EachLv_IdxTable[lv] = NULL;
      }
   }


// 2-2-2. load LBIdx list (sorted by GID)
   H5_SetID_LBIdx = H5Dopen( H5_FileID, "Tree/LBIdx", H5P_DEFAULT );

   if ( H5_SetID_LBIdx < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/LBIdx" );

   H5_Status = H5Dread( H5_SetID_LBIdx, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, LBIdxList_AllLv );
   H5_Status = H5Dclose( H5_SetID_LBIdx );


   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
//    2-2-3. sort the LBIdx list at each level and set the load-balance cut points
#     if ( LOAD_BALANCE != HILBERT )
      if ( NLvRescale != 1 )
      Aux_Message( stderr, "WARNING : please make sure that the patch LBIdx doesn't change when NLvRescale != 1 !!\n" );
#     endif

//    Mis_Heapsort must be called before LB_SetCutPoint in order to get LBIdxList_EachLv_IdxTable
//    (since LB_SetCutPoint will sort LBIdxList_EachLv as well)
//    --> Actually it's not necessary anymore since we now send LBIdx0_AllRank instead of LBIdxList_EachLv into LB_SetCutPoint()
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_EachLv[lv], LBIdxList_EachLv_IdxTable[lv] );

//    prepare LBIdx and load-balance weighting of each **patch group** for LB_SetCutPoint()
      const bool InputLBIdx0AndLoad_Yes = true;
      long   *LBIdx0_AllRank = NULL;
      double *Load_AllRank   = NULL;

      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = new long   [ NPatchTotal[lv] / 8 ];
         Load_AllRank   = new double [ NPatchTotal[lv] / 8 ];

         for (int t=0; t<NPatchTotal[lv]/8; t++)
         {
            LBIdx0_AllRank[t]  = LBIdxList_EachLv[lv][t*8];
            LBIdx0_AllRank[t] -= LBIdx0_AllRank[t] % 8;     // assuming LBIdxList_EachLv is NOT sorted yet
            Load_AllRank  [t]  = 8.0;                       // assuming all patches have the same weighting == 1.0
         }
      }

//    do NOT consider load-balance weighting of particles since at this point we don't have that information
      const double ParWeight_Zero = 0.0;
      LB_SetCutPoint( lv, NPatchTotal[lv]/8, amr->LB->CutPoint[lv], InputLBIdx0AndLoad_Yes, LBIdx0_AllRank, Load_AllRank,
                      ParWeight_Zero );

//    free memory
      if ( MPI_Rank == 0 )
      {
         delete [] LBIdx0_AllRank;
         delete [] Load_AllRank;
      }


//    2-2-4. get the target LBIdx range of each rank
      LoadIdx_Start[lv] = -1;    // -1 --> indicate that it has not been set
      LoadIdx_Stop [lv] = -1;    // must be initialized as <= LoadIdx_Start[lv]

//    skip levels with no patch (otherwise LoadIdx_Start and LoadIdx_Stop will be set incorrectly)
      if ( NPatchTotal[lv] == 0 )   continue;

      for (int t=0; t<NPatchTotal[lv]; t++)
      {
//       set LoadIdx_Start to "the first patch belonging to this rank"
         if (  LoadIdx_Start[lv] == -1  &&  LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) == MPI_Rank  )
            LoadIdx_Start[lv] = t;

//       set LoadIdx_Stop to "the last patch belonging to this rank + 1"
         if ( LoadIdx_Start[lv] != -1 )
         {
            if (  LB_Index2Rank( lv, LBIdxList_EachLv[lv][t], CHECK_ON ) > MPI_Rank  )
            {
               LoadIdx_Stop[lv] = t;
               break;
            }

//          rank owning the last patch needs to be treated separately
            else if ( t == NPatchTotal[lv] - 1 )
            {
               LoadIdx_Stop[lv] = NPatchTotal[lv];
               break;
            }
         }
      }

#     ifdef DEBUG_HDF5
      if ( LoadIdx_Start[lv]%8 != 0  &&  LoadIdx_Start[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Start[%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Start[lv]%8 );

      if ( LoadIdx_Stop [lv]%8 != 0  &&  LoadIdx_Stop[lv] != -1 )
         Aux_Error( ERROR_INFO, "LoadIdx_Stop [%d] = %d --> %%8 != 0 !!\n", lv, LoadIdx_Stop [lv]%8 );

      if (  ( LoadIdx_Start[lv] == -1 && LoadIdx_Stop[lv] != -1 )  ||
            ( LoadIdx_Start[lv] != -1 && LoadIdx_Stop[lv] == -1 )   )
         Aux_Error( ERROR_INFO, "LoadIdx_Start/Stop[%d] = %d/%d !!\n", lv, LoadIdx_Start[lv], LoadIdx_Stop[lv] );
#     endif
   } // for (int lv=0; lv<KeyInfo.NLevel; lv++)

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading load-balance index table ... done\n" );


#  else // #ifdef LOAD_BALANCE


// 2-3. son
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading son table ...\n" );

// allocate memory
   int *SonList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_Son = H5Dopen( H5_FileID, "Tree/Son", H5P_DEFAULT );
   if ( H5_SetID_Son < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Son" );
   H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Son );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading son table ... done\n" );
#  endif // #ifdef LOAD_BALANCE ... else ...


// 2-4. number of particles in each patch
#  ifdef PARTICLE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading particle counts ...\n" );

   hid_t H5_SetID_NPar;

// allocate memory
   int *NParList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_NPar = H5Dopen( H5_FileID, "Tree/NPar", H5P_DEFAULT );
   if ( H5_SetID_NPar < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/NPar" );
   H5_Status = H5Dread( H5_SetID_NPar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, NParList_AllLv );
   H5_Status = H5Dclose( H5_SetID_NPar );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading particle counts ... done\n" );
#  endif // #ifdef PARTICLE


   H5_Status = H5Fclose( H5_FileID );


// 2-5. initialize particle variables
#  ifdef PARTICLE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing particle repository ...\n" );

// 2-5-1. get the total number of paticles in this rank
   long NParThisRank;

#  ifdef LOAD_BALANCE
   NParThisRank = 0;

   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   for (int t=LoadIdx_Start[lv]; t<LoadIdx_Stop[lv]; t++)
   {
#     ifdef DEBUG_HDF5
      if ( t < 0  ||  t >= NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "incorrect load index (%d) !!\n", t );
#     endif

      NParThisRank += NParList_AllLv[ LBIdxList_EachLv_IdxTable[lv][t] + GID_LvStart[lv] ];
   }

#  ifdef DEBUG_HDF5
   long NParAllRank;
   MPI_Reduce( &NParThisRank, &NParAllRank, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   if ( MPI_Rank == 0  &&  NParAllRank != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles (%ld) != expect (%ld) !!\n",
                 NParAllRank, amr->Par->NPar_Active_AllRank );
#  endif

#  else // #ifdef LOAD_BALANCE

   NParThisRank = amr->Par->NPar_Active_AllRank;
#  endif // #ifdef LOAD_BALANCE ... else ...


// 2-5-2. initialize particle repository
   amr->Par->InitRepo( NParThisRank, MPI_NRank );

// reset the total number of particles to be zero
// --> so particle repository is pre-allocated, but it contains no active particle yet
// --> particles will be added later by calling Par->AddOneParticle()
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;


// 2-5-3. calculate the starting global particle indices (i.e., GParID_Offset) for all patches
   long *GParID_Offset = new long [ NPatchAllLv ];

   GParID_Offset[0] = 0;
   for (int t=1; t<NPatchAllLv; t++)   GParID_Offset[t] = GParID_Offset[t-1] + NParList_AllLv[t-1];

#  ifdef DEBUG_HDF5
   if ( GParID_Offset[ NPatchAllLv-1 ] + NParList_AllLv[ NPatchAllLv-1 ] != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles (%ld) != expect (%ld) !!\n",
                 GParID_Offset[ NPatchAllLv-1 ] + NParList_AllLv[ NPatchAllLv-1 ], amr->Par->NPar_Active_AllRank );
#  endif


// 2-5-4. get the maximum number of particles in one patch and allocate an I/O buffer accordingly
   long MaxNParInOnePatch = 0;
   long *NewParList       = NULL;
   real **ParBuf          = NULL;

   for (int t=0; t<NPatchAllLv; t++)   MaxNParInOnePatch = MAX( MaxNParInOnePatch, NParList_AllLv[t] );

   NewParList = new long [MaxNParInOnePatch];

// be careful about using ParBuf returned from Aux_AllocateArray2D, which is set to NULL if MaxNParInOnePatch == 0
// --> for example, accessing ParBuf[0...PAR_NATT_STORED-1] will be illegal when MaxNParInOnePatch == 0
   Aux_AllocateArray2D( ParBuf, PAR_NATT_STORED, MaxNParInOnePatch );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing particle repository ... done\n" );
#  endif // #ifdef PARTICLE



// 3. load and allocate patches (load particles as well if PARTICLE is on)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches and particles ...\n" );

#  ifdef LOAD_BALANCE
   const bool Recursive_No  = false;
#  else
   const bool Recursive_Yes = true;
   int  TRange_Min[3], TRange_Max[3];
#  endif

   char (*FieldName)[MAX_STRING] = new char [NCOMP_TOTAL][MAX_STRING];
   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[NCOMP_TOTAL], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_GridData;

#  ifdef MHD
   char (*FCMagName)[MAX_STRING] = new char [NCOMP_MAG][MAX_STRING];
   hsize_t H5_SetDims_FCMag[4], H5_MemDims_FCMag[4];
   hid_t   H5_SetID_FCMag[NCOMP_MAG], H5_MemID_FCMag[NCOMP_MAG], H5_SpaceID_FCMag[NCOMP_MAG];
#  else
   hid_t *H5_SetID_FCMag   = NULL;
   hid_t *H5_MemID_FCMag   = NULL;
   hid_t *H5_SpaceID_FCMag = NULL;
#  endif // #ifdef MHD ... else ...

#  ifdef PARTICLE
   char (*ParAttName)[MAX_STRING] = new char [PAR_NATT_STORED][MAX_STRING];
   hsize_t H5_SetDims_ParData[1];
   hid_t   H5_SetID_ParData[PAR_NATT_STORED], H5_SpaceID_ParData, H5_GroupID_Particle;
#  else
// define useless variables when PARTICLE is off
   int   *NParList_AllLv     = NULL;
   real **ParBuf             = NULL;
   long  *NewParList         = NULL;
   long  *GParID_Offset      = NULL;
   hid_t *H5_SetID_ParData   = NULL;
   hid_t  H5_SpaceID_ParData = NULL_INT;
   long   NParThisRank       = NULL_INT;
#  endif // #ifdef PARTICLE ... else ...


// 3-1. set the names of all grid fields and particle attributes
   for (int v=0; v<NCOMP_TOTAL; v++)      sprintf( FieldName[v], "%s", FieldLabel[v] );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)        sprintf( FCMagName[v], "%s", MagLabel[v] );
#  endif

#  ifdef PARTICLE
// skip the last PAR_NATT_UNSTORED attributes since we do not store them on disk
   for (int v=0; v<PAR_NATT_STORED; v++)  sprintf( ParAttName[v], "%s", ParAttLabel[v] );
#  endif


// 3-2. initialize relevant HDF5 objects
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PS1;
   H5_SetDims_Field[2] = PS1;
   H5_SetDims_Field[3] = PS1;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

   H5_MemDims_Field[0] = 1;
   H5_MemDims_Field[1] = PS1;
   H5_MemDims_Field[2] = PS1;
   H5_MemDims_Field[3] = PS1;

   H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
   if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemID_Field" );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_SetDims_FCMag[0] = NPatchAllLv;
      for (int t=1; t<4; t++)
      H5_SetDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_SpaceID_FCMag[v] = H5Screate_simple( 4, H5_SetDims_FCMag, NULL );
      if ( H5_SpaceID_FCMag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_SpaceID_FCMag", v );

      H5_MemDims_FCMag[0] = 1;
      for (int t=1; t<4; t++)
      H5_MemDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_MemID_FCMag[v] = H5Screate_simple( 4, H5_MemDims_FCMag, NULL );
      if ( H5_MemID_FCMag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_MemID_FCMag", v );
   }
#  endif // #ifdef MHD

#  ifdef PARTICLE
   H5_SetDims_ParData[0] = amr->Par->NPar_Active_AllRank;
   H5_SpaceID_ParData    = H5Screate_simple( 1, H5_SetDims_ParData, NULL );
   if ( H5_SpaceID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_ParData" );
#  endif


// load data with RESTART_LOAD_NRANK ranks at a time
   for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)
   {
      if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )
      {
//       3-3. open the target datasets just once
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
         if ( H5_FileID < 0 )
            Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
         if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
            if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
         }

#        ifdef MHD
         for (int v=0; v<NCOMP_MAG; v++)
         {
            H5_SetID_FCMag[v] = H5Dopen( H5_GroupID_GridData, FCMagName[v], H5P_DEFAULT );
            if ( H5_SetID_FCMag[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FCMagName[v] );
         }
#        endif

#        ifdef PARTICLE
         H5_GroupID_Particle = H5Gopen( H5_FileID, "Particle", H5P_DEFAULT );
         if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Particle" );

         for (int v=0; v<PAR_NATT_STORED; v++)
         {
            H5_SetID_ParData[v] = H5Dopen( H5_GroupID_Particle, ParAttName[v], H5P_DEFAULT );
            if ( H5_SetID_ParData[v] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", ParAttName[v] );
         }
#        endif


//       3-4. begin to load data
//       3-4-1. load-balance data
#        ifdef LOAD_BALANCE
         int GID0;

         for (int lv=0; lv<KeyInfo.NLevel; lv++)
         {
            if ( MPI_Rank == TRanks )
            Aux_Message( stdout, "      Loading ranks %4d -- %4d, lv %2d ... ",
                         TRanks, MIN(TRanks+RESTART_LOAD_NRANK-1, MPI_NRank-1), lv );

//          loop over all target LBIdx
            for (int t=LoadIdx_Start[lv]; t<LoadIdx_Stop[lv]; t+=8)
            {
#              ifdef DEBUG_HDF5
               if ( t < 0  ||  t >= NPatchTotal[lv]  ||  t%8 != 0 )
                 Aux_Error( ERROR_INFO, "incorrect load index (%d) !!\n", t );
#              endif

//             make sure that we load patch from LocalID == 0
               GID0 = LBIdxList_EachLv_IdxTable[lv][t] - LBIdxList_EachLv_IdxTable[lv][t]%8 + GID_LvStart[lv];

               for (int GID=GID0; GID<GID0+8; GID++)
                  LoadOnePatch( H5_FileID, lv, GID, Recursive_No, NULL, CrList_AllLv,
                                H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                                H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                                NParList_AllLv, ParBuf, NewParList, H5_SetID_ParData, H5_SpaceID_ParData,
                                GParID_Offset, NParThisRank );
            }

//          check if LocalID matches corner
#           ifdef DEBUG_HDF5
            const int PatchScale = PS1*amr->scale[lv];
            int LocalID;

            for (int PID=0; PID<amr->num[lv]; PID++)
            {
               LocalID = PID%8;

               for (int d=0; d<3; d++)
               {
                  if (  TABLE_02( LocalID, 'x'+d, 0, 1 ) != (amr->patch[0][lv][PID]->corner[d]/PatchScale)%2  )
                  {
                     Output_Patch( lv, PID, 0, 0, 0, NULL );
                     Aux_Error( ERROR_INFO, "lv %d, PID %d, LocalID does not match patch corner (patch has been output) !!\n",
                                 lv, PID );
                  }
               }
            }
#           endif

            if ( MPI_Rank == TRanks )  Aux_Message( stdout, "done\n" );
         } // for (int lv=0; lv<KeyInfo.NLevel; lv++)


//       3-4.2. non load-balance data
#        else // #ifdef LOAD_BALANCE

//       set the target range of each rank
         for (int d=0; d<3; d++)
         {
            TRange_Min[d] = MPI_Rank_X[d]*NX0[d]*amr->scale[0];
            TRange_Max[d] = TRange_Min[d] + NX0[d]*amr->scale[0];
         }

//       loop over the corners of all root-level patches (rescale in advance if KeyInfo.NLevel != NLEVEL)
         const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

         for (int GID=0; GID<NPatchTotal[0]; GID++)
         {
            if ( GID % TenPercent == 0 )
            Aux_Message( stdout, "      Rank %2d: %5.1lf%% completed ...\n",
                         MPI_Rank, 100.0*GID/NPatchTotal[0] );

//          load the info and data of the entire patch family recursively if the root patch is within the target range
            if (  CrList_AllLv[GID][0] >= TRange_Min[0]  &&  CrList_AllLv[GID][0] < TRange_Max[0]  &&
                  CrList_AllLv[GID][1] >= TRange_Min[1]  &&  CrList_AllLv[GID][1] < TRange_Max[1]  &&
                  CrList_AllLv[GID][2] >= TRange_Min[2]  &&  CrList_AllLv[GID][2] < TRange_Max[2]     )
               LoadOnePatch( H5_FileID, 0, GID, Recursive_Yes, SonList_AllLv, CrList_AllLv,
                             H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                             H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                             NParList_AllLv, ParBuf, NewParList, H5_SetID_ParData, H5_SpaceID_ParData,
                             GParID_Offset, NParThisRank );
         } // for (int GID=0; GID<NPatchTotal[0]; GID++)

#        endif // #ifdef LOAD_BALANCE ... else ...

//       free resource
         for (int v=0; v<NCOMP_TOTAL; v++)      H5_Status = H5Dclose( H5_SetID_Field[v] );
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG;   v++)      H5_Status = H5Dclose( H5_SetID_FCMag[v] );
#        endif
         H5_Status = H5Gclose( H5_GroupID_GridData );

#        ifdef PARTICLE
         for (int v=0; v<PAR_NATT_STORED; v++)  H5_Status = H5Dclose( H5_SetID_ParData[v] );
         H5_Status = H5Gclose( H5_GroupID_Particle );
#        endif

         H5_Status = H5Fclose( H5_FileID );
      } // if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)

// free HDF5 objects
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_Status = H5Sclose( H5_SpaceID_FCMag[v] );
      H5_Status = H5Sclose( H5_MemID_FCMag  [v] );
   }
#  endif
#  ifdef PARTICLE
   H5_Status = H5Sclose( H5_SpaceID_ParData );
#  endif


// 3-5. record the number of real patches (and LB_IdxList_Real)
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = amr->num[lv];

#     ifdef LOAD_BALANCE
      if ( amr->LB->IdxList_Real         [lv] != NULL )   delete [] amr->LB->IdxList_Real         [lv];
      if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )   delete [] amr->LB->IdxList_Real_IdxTable[lv];

      amr->LB->IdxList_Real         [lv] = new long [ amr->NPatchComma[lv][1] ];
      amr->LB->IdxList_Real_IdxTable[lv] = new int  [ amr->NPatchComma[lv][1] ];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

      Mis_Heapsort( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );
#     endif

//    get the total number of real patches at all ranks
      Mis_GetTotalPatchNumber( lv );
   }


// 3-6. verify that all patches and particles are loaded
#  ifdef DEBUG_HDF5
   int NLoadPatch[KeyInfo.NLevel];
   MPI_Reduce( amr->num, NLoadPatch, KeyInfo.NLevel, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   for (int lv=0; lv<KeyInfo.NLevel; lv++)
   {
      if ( NLoadPatch[lv] != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "Lv %d: total number of loaded patches (%d) != expectation (%d) !!\n",
                       lv, NLoadPatch[lv], NPatchTotal[lv] );
   }

#  ifdef PARTICLE
   if ( amr->Par->NPar_AcPlusInac != NParThisRank )
      Aux_Error( ERROR_INFO, "total number of particles in the repository (%ld) != expect (%ld) !!\n",
                 amr->Par->NPar_AcPlusInac, NParThisRank );
#  endif
#  endif // #ifdfe DEBUG_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading patches and particles ... done\n" );



// 4. close all HDF5 objects and free memory (especially for the variable-length string)
   free ( KeyInfo.CodeVersion );
   free ( KeyInfo.DumpWallTime );

   delete [] FieldName;
   delete [] CrList_AllLv;
#  ifdef MHD
   delete [] FCMagName;
#  endif
#  ifdef LOAD_BALANCE
   delete [] LBIdxList_AllLv;
   for (int lv=0; lv<NLEVEL; lv++)  delete [] LBIdxList_EachLv_IdxTable[lv];
#  else
   delete [] SonList_AllLv;
#  endif
#  ifdef PARTICLE
   delete [] ParAttName;
   delete [] NParList_AllLv;
   delete [] GParID_Offset;
   delete [] NewParList;
   Aux_DeallocateArray2D( ParBuf );
#  endif



// 5-1. improve load balance
// ===================================================================================================================
#  ifdef LOAD_BALANCE

// no need to redistribute all patches again since we already did that when loading data from disks
// --> but note that we have not considerer load-balance weighting of particles yet

// we don't have enough information to calculate the load-balance weighting of particles when
// calling LB_Init_LoadBalance() for the first time
// --> for example, LB_EstimateWorkload_AllPatchGroup()->Par_CollectParticle2OneLevel()->Par_LB_CollectParticle2OneLevel()
//     needs amr->LB->IdxList_Real[], which will be constructed only AFTER calling LB_Init_LoadBalance()
// --> must disable particle weighting (by setting ParWeight==0.0) first

// must not reset load-balance variables (i.e., must adopt ResetLB_No) when calling LB_Init_LoadBalance() for the first time
// since we MUST NOT overwrite IdxList_Real[] and IdxList_Real_IdxList[] already set above
   const double ParWeight_Zero   = 0.0;
   const bool   Redistribute_Yes = true;
   const bool   Redistribute_No  = false;
   const bool   ResetLB_Yes      = true;
   const bool   ResetLB_No       = false;
   const int    AllLv            = -1;

   LB_Init_LoadBalance( Redistribute_No, ParWeight_Zero, ResetLB_No, AllLv );

// redistribute patches again if we want to take into account the load-balance weighting of particles
#  ifdef PARTICLE
   if ( amr->LB->Par_Weight > 0.0 )
   LB_Init_LoadBalance( Redistribute_Yes, amr->LB->Par_Weight, ResetLB_Yes, AllLv );
#  endif



// 5-2. complete all levels for the case without LOAD_BALANCE
// ===================================================================================================================
#  else // #ifdef LOAD_BALANCE

   for (int lv=0; lv<NLEVEL; lv++)
   {
//    construct the relation "father <-> son"
      if ( lv > 0 )     FindFather( lv, 1 );

//    allocate the buffer patches
      Buf_AllocateBufferPatch( amr, lv );

//    set up the BaseP List
      if ( lv == 0 )    Init_RecordBasePatch();

//    set up the BounP_IDMap
      Buf_RecordBoundaryPatch( lv );

//    construct the sibling relation
      SiblingSearch( lv );

//    get the IDs of patches for sending and receiving data between neighbor ranks
      Buf_RecordExchangeDataPatchID( lv );

//    allocate flux arrays on level "lv-1"
      if ( lv > 0  &&  amr->WithFlux )    Flu_AllocateFluxArray( lv-1 );

//    allocate electric arrays on level "lv-1"
#     ifdef MHD
      if ( lv > 0  &&  amr->WithElectric )   MHD_AllocateElectricArray( lv-1 );
#     endif

//    get data for all buffer patches
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );

   } // for (int lv=0; lv<NLEVEL; lv++)

#  endif // #ifdef LOAD_BALANCE ... else ...


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByRestart_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadField
// Description :  Load a single field from the input compound dataset
//
// Note        :  1. This function works for arbitary datatype (int, float, char, 1D array ...)
//                2. Memory must be allocated for FieldPtr in advance with sufficent size (except for "char *")
//                3. For loading a string, which has (type(FieldPtr) = (char *)), the memory must be freed
//                   manually by calling free()
//                4. It can also compare the loaded variables (FieldPtr) with the reference values (ComprPtr)
//                   (perform comparison only if "NCompr > 0")
//                   --> Please make sure that "FieldPtr" and "ComprPtr" point to the same type since we
//                       use the type of "ComprPtr" to typecast "FieldPtr"
//
// Parameter   :  FieldName        : Name of the target field
//                FieldPtr         : Pointer to store the retrieved data
//                H5_SetID_Target  : HDF5 dataset  ID of the target compound variable
//                H5_TypeID_Target : HDF5 datatype ID of the target compound variable
//                Fatal_Nonexist   : Whether or not the nonexistence of the target field is fatal
//                                   --> true  : terminate the program     if the target field cannot be found
//                                       false : display a warning message if the target field cannot be found
//                ComprPtr         : Pointer to store the reference values for comparison
//                NCompr           : Number of elements to be compared
//                Fatal_Compr      : Whether or not the comparison result is fatal
//                                   --> true  : terminate the program     if "FieldPtr[X] != ComprPtr[X]"
//                                       false : display a warning message if "FieldPtr[X] != ComprPtr[X]"
//
// Return      :  Success/fail <-> 0/<0
//-------------------------------------------------------------------------------------------------------
template <typename T>
herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                  const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                  const T *ComprPtr, const int NCompr, const bool Fatal_Compr )
{

// nothing to do if NCompr == 0 (note that in certain circumstances some variables can have zero size)
   if ( NCompr == 0 )   return 0;


#  ifdef DEBUG_HDF5
   if ( NCompr > 0  &&  ComprPtr == NULL )
      Aux_Error( ERROR_INFO, "ComprPtr == NULL for NCompr = %d > 0 !!\n", NCompr );
#  endif


   bool   Check_Pass = true;
   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;


// load
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx >= 0 )
   {
      H5_TypeID_Field  = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize     = H5Tget_size( H5_TypeID_Field );

      H5_TypeID_Load   = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status        = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

      H5_Status        = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldName );

      H5_Status        = H5Tclose( H5_TypeID_Field );
      H5_Status        = H5Tclose( H5_TypeID_Load  );
   } // if ( H5_FieldIdx >= 0 )

   else
   {
      if ( Fatal_Nonexist )
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the restart file !!\n", FieldName );

      else if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist in the restart file !!\n", FieldName );

      return -1;
   } // if ( H5_FieldIdx >= 0 ) ... else ...


// comparison
   char ArrayIdx[10];

// compare strings
   if ( NCompr > 0  &&  typeid(T) == typeid(char) )
   {
      if (  strcmp( *(char**)FieldPtr, (char*)ComprPtr ) != 0  )
      {
         if ( Fatal_Compr )
         {
            Aux_Error( ERROR_INFO, "\"%s\" : RESTART file (%s) != runtime (%s) !!\n",
                       FieldName, *(char**)FieldPtr, (char*)ComprPtr );
            return -2;
         }

         else
         {
            if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : \"%s\" : RESTART file (%s) != runtime (%s) !!\n",
                         FieldName, *(char**)FieldPtr, (char*)ComprPtr );
            Check_Pass = false;
         }
      }
   }

// compare other data types
   else
   for (int t=0; t<NCompr; t++)
   {
      if ( NCompr == 1 )   sprintf( ArrayIdx, "%s", "" );
      else                 sprintf( ArrayIdx, "[%d]", t );

      if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
      {
         if (  typeid(T) == typeid( int)  ||  typeid(T) == typeid( long)  ||
               typeid(T) == typeid(uint)  ||  typeid(T) == typeid(ulong)    )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n",
                          FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               return -2;
            }

            else
            {
               if ( MPI_Rank == 0 )
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n",
                            FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else if (  typeid(T) == typeid(float)  ||  typeid(T) == typeid(double)  )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                          FieldName, ArrayIdx,  ((T*)FieldPtr)[t], ComprPtr[t] );
               return -2;
            }

            else
            {
               if ( MPI_Rank == 0 )
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                            FieldName, ArrayIdx, ((T*)FieldPtr)[t], ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else
         {
            if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : \"%s%s\" : unsupported data type !!\n",
                         FieldName, ArrayIdx );
            Check_Pass = false;
         }

      } // if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
   } // for (int t=0; t<NCompr; t++)

   return (Check_Pass) ? 0 : -3;

} // FUNCTION : LoadField



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatch
// Description :  Allocate and load all fields (and particles if PARTICLE is on) for one patch
//
// Note        :  1. Both leaf and non-leaf patches will store data in the HDF5 output
//                   --> But only leaf patches have particles associated with them
//                2. If "Recursive == true", this function will be invoked recursively to find all children
//                   (and children's children, ...) patches
//
// Parameter   :  H5_FileID          : HDF5 file ID of the restart file
//                lv                 : Target level
//                GID                : Target GID
//                Recursive          : Find all children (and childrens' children, ...) recuresively
//                SonList            : List of son indices
//                                     --> Set only when LOAD_BALANCE is not defined
//                CrList             : List of patch corners
//                H5_SetID_Field     : HDF5 dataset ID for cell-centered grid data
//                H5_SpaceID_Field   : HDF5 dataset dataspace ID for cell-centered grid data
//                H5_MemID_Field     : HDF5 memory dataspace ID for cell-centered grid data
//                H5_SetID_FCMag     : HDF5 dataset ID for face-centered magnetic field
//                H5_SpaceID_FCMag   : HDF5 dataset dataspace ID for face-centered magnetic field
//                H5_MemID_FCMag     : HDF5 memory dataspace ID for face-centered magnetic field
//                NParList           : List of particle counts
//                ParBuf             : I/O buffer for loading particle data from the disk
//                                     --> It must be preallocated with a size equal to the maximum number of
//                                         particles in one patch times the number of particles attributes
//                                         stored on disk
//                                     --> Be careful about using ParBuf, which is set to NULL if it has no elements
//                                         (because of the current implementation of Aux_AllocateArray2D)
//                                         --> For example, accessing ParBuf[0...PAR_NATT_STORED-1] will be illegal when there
//                                             are no particles
//                NewParList         : Array to store the new particle indices
//                                     --> It must be preallocated with a size equal to the maximum number of
//                                         particles in one patch
//                H5_SetID_ParData   : HDF5 dataset ID for particle data
//                H5_SpaceID_ParData : HDF5 dataset dataspace ID for particle data
//                GParID_Offset      : Starting global particle indices for all patches
//                NParThisRank       : Total number of particles in this rank (for check only)
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive,
                   const int *SonList, const int (*CrList)[3],
                   const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                   const hid_t *H5_SetID_FCMag, const hid_t *H5_SpaceID_FCMag, const hid_t *H5_MemID_FCMag,
                   const int *NParList, real **ParBuf, long *NewParList, const hid_t *H5_SetID_ParData,
                   const hid_t H5_SpaceID_ParData, const long *GParID_Offset, const long NParThisRank )
{

   const bool WithData_Yes = true;

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
#  ifdef MHD
   hsize_t H5_Count_FCMag[4], H5_Offset_FCMag[4];
#  endif
   herr_t  H5_Status;
   int     SonGID0, PID;

// allocate patch
   amr->pnew( lv, CrList[GID][0], CrList[GID][1], CrList[GID][2], -1, WithData_Yes, WithData_Yes, WithData_Yes );

   PID = amr->num[lv] - 1;


// determine the subset of dataspace for grid data
   H5_Offset_Field[0] = GID;
   H5_Offset_Field[1] = 0;
   H5_Offset_Field[2] = 0;
   H5_Offset_Field[3] = 0;

   H5_Count_Field [0] = 1;
   H5_Count_Field [1] = PS1;
   H5_Count_Field [2] = PS1;
   H5_Count_Field [3] = PS1;

   H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
   if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );


// load cell-centered intrinsic variables from disk
// --> excluding all derived variables such as gravitational potential and cell-centered B field
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                           amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   }


// load face-centered magnetic field from disk
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
//    determine the subset of dataspace for grid data
      H5_Offset_FCMag[0] = GID;
      H5_Offset_FCMag[1] = 0;
      H5_Offset_FCMag[2] = 0;
      H5_Offset_FCMag[3] = 0;

      H5_Count_FCMag [0] = 1;
      for (int t=1; t<4; t++)
      H5_Count_FCMag [t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_FCMag[v], H5S_SELECT_SET, H5_Offset_FCMag, NULL, H5_Count_FCMag, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the magnetic field %d !!\n", v );

//    load data
      H5_Status = H5Dread( H5_SetID_FCMag[v], H5T_GAMER_REAL, H5_MemID_FCMag[v], H5_SpaceID_FCMag[v], H5P_DEFAULT,
                           amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load magnetic field (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   } // for (int v=0; v<NCOMP_MAG; v++)
#  endif // #ifdef MHD


// load particle data
#  ifdef PARTICLE
   const int NParThisPatch = NParList[GID];

   hsize_t H5_Offset_ParData[1], H5_Count_ParData[1], H5_MemDims_ParData[1];
   hid_t   H5_MemID_ParData;
   real    NewParAtt[PAR_NATT_TOTAL];

   if ( NParThisPatch > 0 )
   {
//    check: only leaf patches can have particles (but note that SonList is NOT set for LOAD_BALANCE)
#     ifndef LOAD_BALANCE
      if ( SonList[GID] != -1 )
         Aux_Error( ERROR_INFO, "particles in a non-leaf patch (lv %d, PID %d, GID %d, SonPID %d, NPar %d) !!\n",
                    lv, PID, GID, SonList[GID], NParThisPatch );
#     endif

//    determine the memory space and the subset of dataspace for particle data
      H5_Offset_ParData [0] = GParID_Offset[GID];
      H5_Count_ParData  [0] = NParThisPatch;
      H5_MemDims_ParData[0] = NParThisPatch;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_ParData, H5S_SELECT_SET, H5_Offset_ParData, NULL, H5_Count_ParData, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the particle data !!\n" );

      H5_MemID_ParData = H5Screate_simple( 1, H5_MemDims_ParData, NULL );
      if ( H5_MemID_ParData < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_ParData" );

//    load particle data from disk
      for (int v=0; v<PAR_NATT_STORED; v++)
      {
//       using ParBuf[v] here is safe since it's NOT called when NParThisPatch == 0
         H5_Status = H5Dread( H5_SetID_ParData[v], H5T_GAMER_REAL, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT,
                              ParBuf[v] );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load a particle attribute (lv %d, GID %d, v %d) !!\n", lv, GID, v );
      }

//    store particles to the particle repository (one particle at a time)
      NewParAtt[PAR_TIME] = Time[0];   // all particles are assumed to be synchronized with the base level

      for (int p=0; p<NParThisPatch; p++)
      {
//       skip the last PAR_NATT_UNSTORED attributes since we do not store them on disk
         for (int v=0; v<PAR_NATT_STORED; v++)  NewParAtt[v] = ParBuf[v][p];

         NewParList[p] = amr->Par->AddOneParticle( NewParAtt );

//       check
         if ( NewParList[p] >= NParThisRank )
            Aux_Error( ERROR_INFO, "New particle ID (%ld) >= maximum allowed value (%ld) !!\n",
                       NewParList[p], NParThisRank );
      } // for (int p=0; p<NParThisPatch )
      
//    link particles to this patch
      const real *PType = amr->Par->Type;
#     ifdef DEBUG_PARTICLE
      const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      char Comment[MAX_STRING];
      sprintf( Comment, "%s, lv %d, PID %d, GID %d, NPar %d", __FUNCTION__, lv, PID, GID, NParThisPatch );
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                           PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#     else
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                           PType );
#     endif

//    free resource
      H5_Status = H5Sclose( H5_MemID_ParData );

   } // if ( NParList[GID] > 0 )
#  endif // #ifdef PARTICLE


// enter the next level
   if ( Recursive )
   {
      if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );

      SonGID0 = SonList[GID];

      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field,
                          H5_SetID_FCMag, H5_SpaceID_FCMag, H5_MemID_FCMag,
                          NParList, ParBuf, NewParList, H5_SetID_ParData, H5_SpaceID_ParData,
                          GParID_Offset, NParThisRank );
      }
   }

} // FUNCTION : LoadOnePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Makefile
// Description :  Load and compare the Makefile_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName      : Restart file name
//                FormatVersion : File format ID
//-------------------------------------------------------------------------------------------------------
void Check_Makefile( const char *FileName, const int FormatVersion )
{

   const bool    Fatal = true;
   const bool NonFatal = false;

   herr_t Status;


// 1. fill in the runtime parameters
   Makefile_t RT;    // RT = RunTime
   FillIn_Makefile( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the Makefile dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/Makefile", H5P_DEFAULT );  // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/Makefile" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   Makefile_t RS;    // RS = ReStart

   LoadField( "Model",                  &RS.Model,                  SID, TID, NonFatal, &RT.Model,                  1,    Fatal );
   LoadField( "Gravity",                &RS.Gravity,                SID, TID, NonFatal, &RT.Gravity,                1,    Fatal );
   LoadField( "Comoving",               &RS.Comoving,               SID, TID, NonFatal, &RT.Comoving,               1,    Fatal );
   LoadField( "Particle",               &RS.Particle,               SID, TID, NonFatal, &RT.Particle,               1,    Fatal );

   LoadField( "UseGPU",                 &RS.UseGPU,                 SID, TID, NonFatal, &RT.UseGPU,                 1, NonFatal );
   LoadField( "GAMER_Debug",            &RS.GAMER_Debug,            SID, TID, NonFatal, &RT.GAMER_Debug,            1, NonFatal );
   LoadField( "BitwiseReproducibility", &RS.BitwiseReproducibility, SID, TID, NonFatal, &RT.BitwiseReproducibility, 1, NonFatal );
   LoadField( "Timing",                 &RS.Timing,                 SID, TID, NonFatal, &RT.Timing,                 1, NonFatal );
   LoadField( "TimingSolver",           &RS.TimingSolver,           SID, TID, NonFatal, &RT.TimingSolver,           1, NonFatal );
   LoadField( "Float8",                 &RS.Float8,                 SID, TID, NonFatal, &RT.Float8,                 1, NonFatal );
   LoadField( "Serial",                 &RS.Serial,                 SID, TID, NonFatal, &RT.Serial,                 1, NonFatal );
   LoadField( "LoadBalance",            &RS.LoadBalance,            SID, TID, NonFatal, &RT.LoadBalance,            1, NonFatal );
   LoadField( "OverlapMPI",             &RS.OverlapMPI,             SID, TID, NonFatal, &RT.OverlapMPI,             1, NonFatal );
   LoadField( "OpenMP",                 &RS.OpenMP,                 SID, TID, NonFatal, &RT.OpenMP,                 1, NonFatal );
   LoadField( "GPU_Arch",               &RS.GPU_Arch,               SID, TID, NonFatal, &RT.GPU_Arch,               1, NonFatal );
   LoadField( "Laohu",                  &RS.Laohu,                  SID, TID, NonFatal, &RT.Laohu,                  1, NonFatal );
   LoadField( "SupportHDF5",            &RS.SupportHDF5,            SID, TID, NonFatal, &RT.SupportHDF5,            1, NonFatal );
   LoadField( "SupportGSL",             &RS.SupportGSL,             SID, TID, NonFatal, &RT.SupportGSL,             1, NonFatal );
   LoadField( "SupportGrackle",         &RS.SupportGrackle,         SID, TID, NonFatal, &RT.SupportGrackle,         1, NonFatal );
   LoadField( "RandomNumber",           &RS.RandomNumber,           SID, TID, NonFatal, &RT.RandomNumber,           1, NonFatal );

   LoadField( "NLevel",                 &RS.NLevel,                 SID, TID, NonFatal, &RT.NLevel,                 1, NonFatal );
   LoadField( "MaxPatch",               &RS.MaxPatch,               SID, TID, NonFatal, &RT.MaxPatch,               1, NonFatal );

#  ifdef GRAVITY
   LoadField( "PotScheme",              &RS.PotScheme,              SID, TID, NonFatal, &RT.PotScheme,              1, NonFatal );
   LoadField( "StorePotGhost",          &RS.StorePotGhost,          SID, TID, NonFatal, &RT.StorePotGhost,          1, NonFatal );
   LoadField( "UnsplitGravity",         &RS.UnsplitGravity,         SID, TID, NonFatal, &RT.UnsplitGravity,         1, NonFatal );
#  endif

#  if   ( MODEL == HYDRO )
   LoadField( "FluScheme",              &RS.FluScheme,              SID, TID, NonFatal, &RT.FluScheme,              1, NonFatal );
#  ifdef LR_SCHEME
   LoadField( "LRScheme",               &RS.LRScheme,               SID, TID, NonFatal, &RT.LRScheme,               1, NonFatal );
#  endif
#  ifdef RSOLVER
   LoadField( "RSolver",                &RS.RSolver,                SID, TID, NonFatal, &RT.RSolver,                1, NonFatal );
#  endif
   LoadField( "DualEnergy",             &RS.DualEnergy,             SID, TID, NonFatal, &RT.DualEnergy,             1, NonFatal );
   LoadField( "Magnetohydrodynamics",   &RS.Magnetohydrodynamics,   SID, TID, NonFatal, &RT.Magnetohydrodynamics,   1,    Fatal );
   LoadField( "CosmicRay",              &RS.CosmicRay,              SID, TID, NonFatal, &RT.CosmicRay,              1,    Fatal );
   LoadField( "EoS",                    &RS.EoS,                    SID, TID, NonFatal, &RT.EoS,                    1, NonFatal );
   LoadField( "BarotropicEoS",          &RS.BarotropicEoS,          SID, TID, NonFatal, &RT.BarotropicEoS,          1, NonFatal );

#  elif ( MODEL == ELBDM )
   LoadField( "ConserveMass",           &RS.ConserveMass,           SID, TID, NonFatal, &RT.ConserveMass,           1, NonFatal );
   LoadField( "Laplacian4th",           &RS.Laplacian4th,           SID, TID, NonFatal, &RT.Laplacian4th,           1, NonFatal );
   LoadField( "SelfInteraction4",       &RS.SelfInteraction4,       SID, TID, NonFatal, &RT.SelfInteraction4,       1, NonFatal );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   LoadField( "StoreParAcc",            &RS.StoreParAcc,            SID, TID, NonFatal, &RT.StoreParAcc,            1, NonFatal );
   LoadField( "StarFormation",          &RS.StarFormation,          SID, TID, NonFatal, &RT.StarFormation,          1, NonFatal );
   if ( FormatVersion >= 2300 )
   LoadField( "Par_NAttUser",           &RS.Par_NAttUser,           SID, TID, NonFatal, &RT.Par_NAttUser,           1,    Fatal );
   else
   LoadField( "Par_NAttUser",           &RS.Par_NAttUser,           SID, TID, NonFatal, &RT.Par_NAttUser,           1, NonFatal );
#  endif


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : Check_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_SymConst
// Description :  Load and compare the SymConst_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName      : Restart file name
//                FormatVersion : File format ID
//-------------------------------------------------------------------------------------------------------
void Check_SymConst( const char *FileName, const int FormatVersion )
{

   const bool    Fatal = true;
   const bool NonFatal = false;

   herr_t Status;


// 1. fill in the runtime parameters
   SymConst_t RT;    // RT = RunTime
   FillIn_SymConst( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the SymConst dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/SymConst", H5P_DEFAULT );  // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/SymConst" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   SymConst_t RS;    // RS = ReStart


   LoadField( "NCompFluid",           &RS.NCompFluid,           SID, TID, NonFatal, &RT.NCompFluid,            1,    Fatal );
   LoadField( "NCompPassive",         &RS.NCompPassive,         SID, TID, NonFatal, &RT.NCompPassive,          1,    Fatal );
   LoadField( "PatchSize",            &RS.PatchSize,            SID, TID, NonFatal, &RT.PatchSize,             1,    Fatal );
   LoadField( "Flu_NIn",              &RS.Flu_NIn,              SID, TID, NonFatal, &RT.Flu_NIn,               1, NonFatal );
   LoadField( "Flu_NOut",             &RS.Flu_NOut,             SID, TID, NonFatal, &RT.Flu_NOut,              1, NonFatal );
   LoadField( "Flu_NIn_T",            &RS.Flu_NIn_T,            SID, TID, NonFatal, &RT.Flu_NIn_T,             1, NonFatal );
   LoadField( "Flu_NIn_S",            &RS.Flu_NIn_S,            SID, TID, NonFatal, &RT.Flu_NIn_S,             1, NonFatal );
   LoadField( "Flu_NOut_S",           &RS.Flu_NOut_S,           SID, TID, NonFatal, &RT.Flu_NOut_S,            1, NonFatal );
   LoadField( "NFluxFluid",           &RS.NFluxFluid,           SID, TID, NonFatal, &RT.NFluxFluid,            1, NonFatal );
   LoadField( "NFluxPassive",         &RS.NFluxPassive,         SID, TID, NonFatal, &RT.NFluxPassive,          1, NonFatal );
   LoadField( "Flu_GhostSize",        &RS.Flu_GhostSize,        SID, TID, NonFatal, &RT.Flu_GhostSize,         1, NonFatal );
   LoadField( "Flu_Nxt",              &RS.Flu_Nxt,              SID, TID, NonFatal, &RT.Flu_Nxt,               1, NonFatal );
   LoadField( "Debug_HDF5",           &RS.Debug_HDF5,           SID, TID, NonFatal, &RT.Debug_HDF5,            1, NonFatal );
   LoadField( "SibOffsetNonperiodic", &RS.SibOffsetNonperiodic, SID, TID, NonFatal, &RT.SibOffsetNonperiodic,  1, NonFatal );
#  ifdef LOAD_BALANCE
   LoadField( "SonOffsetLB",          &RS.SonOffsetLB,          SID, TID, NonFatal, &RT.SonOffsetLB,           1, NonFatal );
#  endif
   LoadField( "TinyNumber",           &RS.TinyNumber,           SID, TID, NonFatal, &RT.TinyNumber,            1, NonFatal );
   LoadField( "HugeNumber",           &RS.HugeNumber,           SID, TID, NonFatal, &RT.HugeNumber,            1, NonFatal );
   LoadField( "MaxError",             &RS.MaxError,             SID, TID, NonFatal, &RT.MaxError,              1, NonFatal );

#  ifdef GRAVITY
   LoadField( "Gra_NIn",              &RS.Gra_NIn,              SID, TID, NonFatal, &RT.Gra_NIn,               1, NonFatal );
   LoadField( "Pot_GhostSize",        &RS.Pot_GhostSize,        SID, TID, NonFatal, &RT.Pot_GhostSize,         1, NonFatal );
   LoadField( "Gra_GhostSize",        &RS.Gra_GhostSize,        SID, TID, NonFatal, &RT.Gra_GhostSize,         1, NonFatal );
   LoadField( "Rho_GhostSize",        &RS.Rho_GhostSize,        SID, TID, NonFatal, &RT.Rho_GhostSize,         1, NonFatal );
   LoadField( "Pot_Nxt",              &RS.Pot_Nxt,              SID, TID, NonFatal, &RT.Pot_Nxt,               1, NonFatal );
   LoadField( "Gra_Nxt",              &RS.Gra_Nxt,              SID, TID, NonFatal, &RT.Gra_Nxt,               1, NonFatal );
   LoadField( "Rho_Nxt",              &RS.Rho_Nxt,              SID, TID, NonFatal, &RT.Rho_Nxt,               1, NonFatal );
#  ifdef UNSPLIT_GRAVITY
   LoadField( "USG_GhostSizeF",       &RS.USG_GhostSizeF,       SID, TID, NonFatal, &RT.USG_GhostSizeF,        1, NonFatal );
   LoadField( "USG_GhostSizeG",       &RS.USG_GhostSizeG,       SID, TID, NonFatal, &RT.USG_GhostSizeG,        1, NonFatal );
   LoadField( "USG_NxtF",             &RS.USG_NxtF,             SID, TID, NonFatal, &RT.USG_NxtF,              1, NonFatal );
   LoadField( "USG_NxtG",             &RS.USG_NxtG,             SID, TID, NonFatal, &RT.USG_NxtG,              1, NonFatal );
#  endif
   LoadField( "ExtPot_BlockSize",     &RS.ExtPot_BlockSize,     SID, TID, NonFatal, &RT.ExtPot_BlockSize,      1, NonFatal );
   LoadField( "Gra_BlockSize",        &RS.Gra_BlockSize,        SID, TID, NonFatal, &RT.Gra_BlockSize,         1, NonFatal );
   LoadField( "ExtPotNAuxMax",        &RS.ExtPotNAuxMax,        SID, TID, NonFatal, &RT.ExtPotNAuxMax,         1, NonFatal );
   LoadField( "ExtAccNAuxMax",        &RS.ExtAccNAuxMax,        SID, TID, NonFatal, &RT.ExtAccNAuxMax,         1, NonFatal );
   LoadField( "ExtPotNGeneMax",       &RS.ExtPotNGeneMax,       SID, TID, NonFatal, &RT.ExtPotNGeneMax,        1, NonFatal );
#  if   ( POT_SCHEME == SOR )
   LoadField( "Pot_BlockSize_z",      &RS.Pot_BlockSize_z,      SID, TID, NonFatal, &RT.Pot_BlockSize_z,       1, NonFatal );
   LoadField( "UsePSolver_10to14",    &RS.UsePSolver_10to14,    SID, TID, NonFatal, &RT.UsePSolver_10to14,     1, NonFatal );
   LoadField( "SOR_RhoShared",        &RS.SOR_RhoShared,        SID, TID, NonFatal, &RT.SOR_RhoShared,         1, NonFatal );
   LoadField( "SOR_CPotShared",       &RS.SOR_CPotShared,       SID, TID, NonFatal, &RT.SOR_CPotShared,        1, NonFatal );
   LoadField( "SOR_UseShuffle",       &RS.SOR_UseShuffle,       SID, TID, NonFatal, &RT.SOR_UseShuffle,        1, NonFatal );
   LoadField( "SOR_UsePadding",       &RS.SOR_UsePadding,       SID, TID, NonFatal, &RT.SOR_UsePadding,        1, NonFatal );
   LoadField( "SOR_ModReduction",     &RS.SOR_ModReduction,     SID, TID, NonFatal, &RT.SOR_ModReduction,      1, NonFatal );
#  elif ( POT_SCHEME == MG )
   LoadField( "Pot_BlockSize_x",      &RS.Pot_BlockSize_x,      SID, TID, NonFatal, &RT.Pot_BlockSize_x,       1, NonFatal );
#  endif
#  endif // #ifdef GRAVITY

#  ifdef PARTICLE
   if ( FormatVersion >= 2300 )
   LoadField( "Par_NAttStored",       &RS.Par_NAttStored,       SID, TID, NonFatal, &RT.Par_NAttStored,        1,    Fatal );
   else
   LoadField( "Par_NAttStored",       &RS.Par_NAttStored,       SID, TID, NonFatal, &RT.Par_NAttStored,        1, NonFatal );
   LoadField( "RhoExt_GhostSize",     &RS.RhoExt_GhostSize,     SID, TID, NonFatal, &RT.RhoExt_GhostSize,      1, NonFatal );
   LoadField( "Debug_Particle",       &RS.Debug_Particle,       SID, TID, NonFatal, &RT.Debug_Particle,        1, NonFatal );
   LoadField( "ParList_GrowthFactor", &RS.ParList_GrowthFactor, SID, TID, NonFatal, &RT.ParList_GrowthFactor,  1, NonFatal );
   LoadField( "ParList_ReduceFactor", &RS.ParList_ReduceFactor, SID, TID, NonFatal, &RT.ParList_ReduceFactor,  1, NonFatal );
#  endif

   LoadField( "BitRep_Flux",          &RS.BitRep_Flux,          SID, TID, NonFatal, &RT.BitRep_Flux,           1, NonFatal );
#  ifdef MHD
   LoadField( "BitRep_Electric",      &RS.BitRep_Electric,      SID, TID, NonFatal, &RT.BitRep_Electric,       1, NonFatal );
#  endif

#  if   ( MODEL == HYDRO )
   LoadField( "Flu_BlockSize_x",      &RS.Flu_BlockSize_x,      SID, TID, NonFatal, &RT.Flu_BlockSize_x,       1, NonFatal );
   LoadField( "Flu_BlockSize_y",      &RS.Flu_BlockSize_y,      SID, TID, NonFatal, &RT.Flu_BlockSize_y,       1, NonFatal );
   LoadField( "CheckNegativeInFluid", &RS.CheckNegativeInFluid, SID, TID, NonFatal, &RT.CheckNegativeInFluid,  1, NonFatal );
   LoadField( "CharReconstruction",   &RS.CharReconstruction,   SID, TID, NonFatal, &RT.CharReconstruction,    1, NonFatal );
   LoadField( "LR_Eint",              &RS.LR_Eint,              SID, TID, NonFatal, &RT.LR_Eint,               1, NonFatal );
   LoadField( "CheckIntermediate",    &RS.CheckIntermediate,    SID, TID, NonFatal, &RT.CheckIntermediate,     1, NonFatal );
   LoadField( "HLL_NoRefState",       &RS.HLL_NoRefState,       SID, TID, NonFatal, &RT.HLL_NoRefState,        1, NonFatal );
   LoadField( "HLL_IncludeAllWaves",  &RS.HLL_IncludeAllWaves,  SID, TID, NonFatal, &RT.HLL_IncludeAllWaves,   1, NonFatal );
   LoadField( "HLLC_WaveSpeed",       &RS.HLLC_WaveSpeed,       SID, TID, NonFatal, &RT.HLLC_WaveSpeed,        1, NonFatal );
   LoadField( "HLLE_WaveSpeed",       &RS.HLLE_WaveSpeed,       SID, TID, NonFatal, &RT.HLLE_WaveSpeed,        1, NonFatal );
#  ifdef MHD
   LoadField( "HLLD_WaveSpeed",       &RS.HLLD_WaveSpeed,       SID, TID, NonFatal, &RT.HLLD_WaveSpeed,        1, NonFatal );
#  endif
#  ifdef N_FC_VAR
   LoadField( "N_FC_Var",             &RS.N_FC_Var,             SID, TID, NonFatal, &RT.N_FC_Var,              1, NonFatal );
#  endif
#  ifdef N_SLOPE_PPM
   LoadField( "N_Slope_PPM",          &RS.N_Slope_PPM,          SID, TID, NonFatal, &RT.N_Slope_PPM,           1, NonFatal );
#  endif
#  ifdef MHD
   LoadField( "EulerY",               &RS.EulerY,               SID, TID, NonFatal, &RT.EulerY,                1, NonFatal );
#  endif
   LoadField( "EoSNAuxMax",           &RS.EoSNAuxMax,           SID, TID, NonFatal, &RT.EoSNAuxMax,            1, NonFatal );
   LoadField( "EoSNTableMax",         &RS.EoSNTableMax,         SID, TID, NonFatal, &RT.EoSNTableMax,          1, NonFatal );

#  elif  ( MODEL == ELBDM )
   LoadField( "Flu_BlockSize_x",      &RS.Flu_BlockSize_x,      SID, TID, NonFatal, &RT.Flu_BlockSize_x,       1, NonFatal );
   LoadField( "Flu_BlockSize_y",      &RS.Flu_BlockSize_y,      SID, TID, NonFatal, &RT.Flu_BlockSize_y,       1, NonFatal );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   LoadField( "dt_Flu_BlockSize",     &RS.dt_Flu_BlockSize,     SID, TID, NonFatal, &RT.dt_Flu_BlockSize,      1, NonFatal );
   LoadField( "dt_Flu_UseShuffle",    &RS.dt_Flu_UseShuffle,    SID, TID, NonFatal, &RT.dt_Flu_UseShuffle,     1, NonFatal );
#  ifdef GRAVITY
   LoadField( "dt_Gra_BlockSize",     &RS.dt_Gra_BlockSize,     SID, TID, NonFatal, &RT.dt_Gra_BlockSize,      1, NonFatal );
   LoadField( "dt_Gra_UseShuffle",    &RS.dt_Gra_UseShuffle,    SID, TID, NonFatal, &RT.dt_Gra_UseShuffle,     1, NonFatal );
#  endif

   LoadField( "Src_BlockSize",        &RS.Src_BlockSize,        SID, TID, NonFatal, &RT.Src_BlockSize,         1, NonFatal );
   LoadField( "Src_GhostSize",        &RS.Src_GhostSize,        SID, TID, NonFatal, &RT.Src_GhostSize,         1, NonFatal );
   LoadField( "Src_Nxt",              &RS.Src_Nxt,              SID, TID, NonFatal, &RT.Src_Nxt,               1, NonFatal );
   LoadField( "Src_NAuxDlep",         &RS.Src_NAuxDlep,         SID, TID, NonFatal, &RT.Src_NAuxDlep,          1, NonFatal );
   LoadField( "Src_DlepProfNVar",     &RS.Src_DlepProfNVar,     SID, TID, NonFatal, &RT.Src_DlepProfNVar,      1, NonFatal );
   LoadField( "Src_DlepProfNBinMax",  &RS.Src_DlepProfNBinMax,  SID, TID, NonFatal, &RT.Src_DlepProfNBinMax,   1, NonFatal );
   LoadField( "Src_NAuxUser",         &RS.Src_NAuxUser,         SID, TID, NonFatal, &RT.Src_NAuxUser,          1, NonFatal );

   LoadField( "Der_GhostSize",        &RS.Der_GhostSize,        SID, TID, NonFatal, &RT.Der_GhostSize,         1, NonFatal );
   LoadField( "Der_Nxt",              &RS.Der_Nxt,              SID, TID, NonFatal, &RT.Der_Nxt,               1, NonFatal );
   LoadField( "Der_NOut_Max",         &RS.Der_NOut_Max,         SID, TID, NonFatal, &RT.Der_NOut_Max,          1, NonFatal );


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : Check_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_InputPara
// Description :  Load and compare the InputPara_t structure (runtime vs. restart file)
//
// Note        :  1. Unnecessary comparisons are avoided
//                2. Program may be terminated if any fatal difference is found
//
// Parameter   :  FileName      : Restart file name
//                FormatVersion : File format ID
//-------------------------------------------------------------------------------------------------------
void Check_InputPara( const char *FileName, const int FormatVersion )
{

   const bool    Fatal = true;
   const bool NonFatal = false;
   const int  N1       = MAX_LEVEL;
   const int  NP       = NCOMP_PASSIVE;
   const int *NullPtr  = NULL;

   herr_t Status;


// 1. fill in the runtime parameters
   InputPara_t RT;   // RT = RunTime
   FillIn_InputPara( RT );


// 2. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 3. load the InputPara dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/InputPara", H5P_DEFAULT ); // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 4. load and compare all target fields one-by-one
   InputPara_t RS;   // RS = ReStart

// simulation scale
   LoadField( "BoxSize",                 &RS.BoxSize,                 SID, TID, NonFatal, &RT.BoxSize,                  1,    Fatal );
   LoadField( "NX0_Tot",                  RS.NX0_Tot,                 SID, TID, NonFatal,  RT.NX0_Tot,                  3,    Fatal );
   LoadField( "MPI_NRank",               &RS.MPI_NRank,               SID, TID, NonFatal, &RT.MPI_NRank,                1, NonFatal );
   LoadField( "MPI_NRank_X",              RS.MPI_NRank_X,             SID, TID, NonFatal,  RT.MPI_NRank_X,              3, NonFatal );
   LoadField( "OMP_NThread",             &RS.OMP_NThread,             SID, TID, NonFatal, &RT.OMP_NThread,              1, NonFatal );
   LoadField( "EndT",                    &RS.EndT,                    SID, TID, NonFatal, &RT.EndT,                     1, NonFatal );
   LoadField( "EndStep",                 &RS.EndStep,                 SID, TID, NonFatal, &RT.EndStep,                  1, NonFatal );

// test problems
   LoadField( "TestProb_ID",             &RS.TestProb_ID,             SID, TID, NonFatal, &RT.TestProb_ID,              1, NonFatal );

// code units
   LoadField( "Opt__Unit",               &RS.Opt__Unit,               SID, TID, NonFatal, &RT.Opt__Unit,                1, NonFatal );
   LoadField( "Unit_L",                  &RS.Unit_L,                  SID, TID, NonFatal, &RT.Unit_L,                   1, NonFatal );
   LoadField( "Unit_M",                  &RS.Unit_M,                  SID, TID, NonFatal, &RT.Unit_M,                   1, NonFatal );
   LoadField( "Unit_T",                  &RS.Unit_T,                  SID, TID, NonFatal, &RT.Unit_T,                   1, NonFatal );
   LoadField( "Unit_V",                  &RS.Unit_V,                  SID, TID, NonFatal, &RT.Unit_V,                   1, NonFatal );
   LoadField( "Unit_D",                  &RS.Unit_D,                  SID, TID, NonFatal, &RT.Unit_D,                   1, NonFatal );
   LoadField( "Unit_E",                  &RS.Unit_E,                  SID, TID, NonFatal, &RT.Unit_E,                   1, NonFatal );
   LoadField( "Unit_P",                  &RS.Unit_P,                  SID, TID, NonFatal, &RT.Unit_P,                   1, NonFatal );
#  ifdef MHD
   LoadField( "Unit_B",                  &RS.Unit_B,                  SID, TID, NonFatal, &RT.Unit_B,                   1, NonFatal );
#  endif

// boundary condition
   LoadField( "Opt__BC_Flu",              RS.Opt__BC_Flu,             SID, TID, NonFatal,  RT.Opt__BC_Flu,              6, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__BC_Pot",             &RS.Opt__BC_Pot,             SID, TID, NonFatal, &RT.Opt__BC_Pot,              1, NonFatal );
   LoadField( "GFunc_Coeff0",            &RS.GFunc_Coeff0,            SID, TID, NonFatal, &RT.GFunc_Coeff0,             1, NonFatal );
#  endif

// particle
#  ifdef PARTICLE
   LoadField( "Par_Init",                &RS.Par_Init,                SID, TID, NonFatal, &RT.Par_Init,                 1, NonFatal );
   LoadField( "Par_ICFormat",            &RS.Par_ICFormat,            SID, TID, NonFatal, &RT.Par_ICFormat,             1, NonFatal );
   LoadField( "Par_ICMass",              &RS.Par_ICMass,              SID, TID, NonFatal, &RT.Par_ICMass,               1, NonFatal );
   LoadField( "Par_Interp",              &RS.Par_Interp,              SID, TID, NonFatal, &RT.Par_Interp,               1, NonFatal );
   LoadField( "Par_Integ",               &RS.Par_Integ,               SID, TID, NonFatal, &RT.Par_Integ,                1, NonFatal );
   LoadField( "Par_ImproveAcc",          &RS.Par_ImproveAcc,          SID, TID, NonFatal, &RT.Par_ImproveAcc,           1, NonFatal );
   LoadField( "Par_PredictPos",          &RS.Par_PredictPos,          SID, TID, NonFatal, &RT.Par_PredictPos,           1, NonFatal );
   LoadField( "Par_RemoveCell",          &RS.Par_RemoveCell,          SID, TID, NonFatal, &RT.Par_RemoveCell,           1, NonFatal );
   LoadField( "Par_GhostSize",           &RS.Par_GhostSize,           SID, TID, NonFatal, &RT.Par_GhostSize,            1, NonFatal );
#  endif

// cosmology
#  ifdef COMOVING
   LoadField( "A_Init",                  &RS.A_Init,                  SID, TID, NonFatal, &RT.A_Init,                   1, NonFatal );
   LoadField( "OmegaM0",                 &RS.OmegaM0,                 SID, TID, NonFatal, &RT.OmegaM0,                  1, NonFatal );
   LoadField( "Hubble0",                 &RS.Hubble0,                 SID, TID, NonFatal, &RT.Hubble0,                  1, NonFatal );
#  endif

// time-step determination
   LoadField( "Dt__Max",                 &RS.Dt__Max,                 SID, TID, NonFatal, &RT.Dt__Max,                  1, NonFatal );
   LoadField( "Dt__Fluid",               &RS.Dt__Fluid,               SID, TID, NonFatal, &RT.Dt__Fluid,                1, NonFatal );
   LoadField( "Dt__FluidInit",           &RS.Dt__FluidInit,           SID, TID, NonFatal, &RT.Dt__FluidInit,            1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Dt__Gravity",             &RS.Dt__Gravity,             SID, TID, NonFatal, &RT.Dt__Gravity,              1, NonFatal );
#  endif
#  if ( MODEL == ELBDM )
   LoadField( "Dt__Phase",               &RS.Dt__Phase,               SID, TID, NonFatal, &RT.Dt__Phase,                1, NonFatal );
#  endif
#  ifdef PARTICLE
   LoadField( "Dt__ParVel",              &RS.Dt__ParVel,              SID, TID, NonFatal, &RT.Dt__ParVel,               1, NonFatal );
   LoadField( "Dt__ParVelMax",           &RS.Dt__ParVelMax,           SID, TID, NonFatal, &RT.Dt__ParVelMax,            1, NonFatal );
   LoadField( "Dt__ParAcc",              &RS.Dt__ParAcc,              SID, TID, NonFatal, &RT.Dt__ParAcc,               1, NonFatal );
#  endif
#  ifdef COMOVING
   LoadField( "Dt__MaxDeltaA",           &RS.Dt__MaxDeltaA,           SID, TID, NonFatal, &RT.Dt__MaxDeltaA,            1, NonFatal );
#  endif
   LoadField( "Dt__SyncParentLv",        &RS.Dt__SyncParentLv,        SID, TID, NonFatal, &RT.Dt__SyncParentLv,         1, NonFatal );
   LoadField( "Dt__SyncChildrenLv",      &RS.Dt__SyncChildrenLv,      SID, TID, NonFatal, &RT.Dt__SyncChildrenLv,       1, NonFatal );
   LoadField( "Opt__DtUser",             &RS.Opt__DtUser,             SID, TID, NonFatal, &RT.Opt__DtUser,              1, NonFatal );
   LoadField( "Opt__DtLevel",            &RS.Opt__DtLevel,            SID, TID, NonFatal, &RT.Opt__DtLevel,             1, NonFatal );
   LoadField( "Opt__RecordDt",           &RS.Opt__RecordDt,           SID, TID, NonFatal, &RT.Opt__RecordDt,            1, NonFatal );
   LoadField( "AutoReduceDt",            &RS.AutoReduceDt,            SID, TID, NonFatal, &RT.AutoReduceDt,             1, NonFatal );
   LoadField( "AutoReduceDtFactor",      &RS.AutoReduceDtFactor,      SID, TID, NonFatal, &RT.AutoReduceDtFactor,       1, NonFatal );
   LoadField( "AutoReduceDtFactorMin",   &RS.AutoReduceDtFactorMin,   SID, TID, NonFatal, &RT.AutoReduceDtFactorMin,    1, NonFatal );


// domain refinement
   LoadField( "RegridCount",             &RS.RegridCount,             SID, TID, NonFatal, &RT.RegridCount,              1, NonFatal );
   LoadField( "FlagBufferSize",          &RS.FlagBufferSize,          SID, TID, NonFatal, &RT.FlagBufferSize,           1, NonFatal );
   LoadField( "FlagBufferSizeMaxM1Lv",   &RS.FlagBufferSizeMaxM1Lv,   SID, TID, NonFatal, &RT.FlagBufferSizeMaxM1Lv,    1, NonFatal );
   LoadField( "FlagBufferSizeMaxM2Lv",   &RS.FlagBufferSizeMaxM2Lv,   SID, TID, NonFatal, &RT.FlagBufferSizeMaxM2Lv,    1, NonFatal );
   LoadField( "MaxLevel",                &RS.MaxLevel,                SID, TID, NonFatal, &RT.MaxLevel,                 1, NonFatal );
   LoadField( "Opt__Flag_Rho",           &RS.Opt__Flag_Rho,           SID, TID, NonFatal, &RT.Opt__Flag_Rho,            1, NonFatal );
   LoadField( "Opt__Flag_RhoGradient",   &RS.Opt__Flag_RhoGradient,   SID, TID, NonFatal, &RT.Opt__Flag_RhoGradient,    1, NonFatal );
#  if ( MODEL == HYDRO )
   LoadField( "Opt__Flag_PresGradient",  &RS.Opt__Flag_PresGradient,  SID, TID, NonFatal, &RT.Opt__Flag_PresGradient,   1, NonFatal );
   LoadField( "Opt__Flag_Vorticity",     &RS.Opt__Flag_Vorticity,     SID, TID, NonFatal, &RT.Opt__Flag_Vorticity,      1, NonFatal );
   LoadField( "Opt__Flag_Jeans",         &RS.Opt__Flag_Jeans,         SID, TID, NonFatal, &RT.Opt__Flag_Jeans,          1, NonFatal );
#  ifdef MHD
   LoadField( "Opt__Flag_Current",       &RS.Opt__Flag_Current,       SID, TID, NonFatal, &RT.Opt__Flag_Current,        1, NonFatal );
#  endif
#  endif
#  if ( MODEL == ELBDM )
   LoadField( "Opt__Flag_EngyDensity",   &RS.Opt__Flag_EngyDensity,   SID, TID, NonFatal, &RT.Opt__Flag_EngyDensity,    1, NonFatal );
#  endif
   LoadField( "Opt__Flag_LohnerDens",    &RS.Opt__Flag_LohnerDens,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerDens,     1, NonFatal );
#  if ( MODEL == HYDRO )
   LoadField( "Opt__Flag_LohnerEngy",    &RS.Opt__Flag_LohnerEngy,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerEngy,     1, NonFatal );
   LoadField( "Opt__Flag_LohnerPres",    &RS.Opt__Flag_LohnerPres,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerPres,     1, NonFatal );
   LoadField( "Opt__Flag_LohnerTemp",    &RS.Opt__Flag_LohnerTemp,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerTemp,     1, NonFatal );
#  endif
   LoadField( "Opt__Flag_LohnerForm",    &RS.Opt__Flag_LohnerForm,    SID, TID, NonFatal, &RT.Opt__Flag_LohnerForm,     1, NonFatal );
   LoadField( "Opt__Flag_User",          &RS.Opt__Flag_User,          SID, TID, NonFatal, &RT.Opt__Flag_User,           1, NonFatal );
   LoadField( "Opt__Flag_User_Num",      &RS.Opt__Flag_User_Num,      SID, TID, NonFatal, &RT.Opt__Flag_User_Num,       1, NonFatal );
   LoadField( "Opt__Flag_Region",        &RS.Opt__Flag_Region,        SID, TID, NonFatal, &RT.Opt__Flag_Region,         1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Flag_NParPatch",     &RS.Opt__Flag_NParPatch,     SID, TID, NonFatal, &RT.Opt__Flag_NParPatch,      1, NonFatal );
   LoadField( "Opt__Flag_NParCell",      &RS.Opt__Flag_NParCell,      SID, TID, NonFatal, &RT.Opt__Flag_NParCell,       1, NonFatal );
   LoadField( "Opt__Flag_ParMassCell",   &RS.Opt__Flag_ParMassCell,   SID, TID, NonFatal, &RT.Opt__Flag_ParMassCell,    1, NonFatal );
#  endif
   LoadField( "Opt__NoFlagNearBoundary", &RS.Opt__NoFlagNearBoundary, SID, TID, NonFatal, &RT.Opt__NoFlagNearBoundary,  1, NonFatal );
   LoadField( "Opt__PatchCount",         &RS.Opt__PatchCount,         SID, TID, NonFatal, &RT.Opt__PatchCount,          1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__ParticleCount",      &RS.Opt__ParticleCount,      SID, TID, NonFatal, &RT.Opt__ParticleCount,       1, NonFatal );
#  endif
   LoadField( "Opt__ReuseMemory",        &RS.Opt__ReuseMemory,        SID, TID, NonFatal, &RT.Opt__ReuseMemory,         1, NonFatal );
   LoadField( "Opt__MemoryPool",         &RS.Opt__MemoryPool,         SID, TID, NonFatal, &RT.Opt__MemoryPool,          1, NonFatal );

// load balance
#  ifdef LOAD_BALANCE
   LoadField( "LB_WLI_Max",              &RS.LB_WLI_Max,              SID, TID, NonFatal, &RT.LB_WLI_Max,               1, NonFatal );
#  ifdef PARTICLE
   LoadField( "LB_Par_Weight",           &RS.LB_Par_Weight,           SID, TID, NonFatal, &RT.LB_Par_Weight,            1, NonFatal );
#  endif
   LoadField( "Opt__RecordLoadBalance",  &RS.Opt__RecordLoadBalance,  SID, TID, NonFatal, &RT.Opt__RecordLoadBalance,   1, NonFatal );
#  endif
   LoadField( "Opt__MinimizeMPIBarrier", &RS.Opt__MinimizeMPIBarrier, SID, TID, NonFatal, &RT.Opt__MinimizeMPIBarrier,  1, NonFatal );

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   LoadField( "Gamma",                   &RS.Gamma,                   SID, TID, NonFatal, &RT.Gamma,                    1, NonFatal );
   LoadField( "MolecularWeight",         &RS.MolecularWeight,         SID, TID, NonFatal, &RT.MolecularWeight,          1, NonFatal );
   LoadField( "IsoTemp",                 &RS.IsoTemp,                 SID, TID, NonFatal, &RT.IsoTemp,                  1, NonFatal );
   LoadField( "MinMod_Coeff",            &RS.MinMod_Coeff,            SID, TID, NonFatal, &RT.MinMod_Coeff,             1, NonFatal );
   LoadField( "Opt__LR_Limiter",         &RS.Opt__LR_Limiter,         SID, TID, NonFatal, &RT.Opt__LR_Limiter,          1, NonFatal );
   LoadField( "Opt__1stFluxCorr",        &RS.Opt__1stFluxCorr,        SID, TID, NonFatal, &RT.Opt__1stFluxCorr,         1, NonFatal );
   LoadField( "Opt__1stFluxCorrScheme",  &RS.Opt__1stFluxCorrScheme,  SID, TID, NonFatal, &RT.Opt__1stFluxCorrScheme,   1, NonFatal );
#  endif

// ELBDM solvers
#  if ( MODEL == ELBDM )
   LoadField( "ELBDM_Mass",              &RS.ELBDM_Mass,              SID, TID, NonFatal, &RT.ELBDM_Mass,               1, NonFatal );
   LoadField( "ELBDM_PlanckConst",       &RS.ELBDM_PlanckConst,       SID, TID, NonFatal, &RT.ELBDM_PlanckConst,        1, NonFatal );
#  ifdef QUARTIC_SELF_INTERACTION
   LoadField( "ELBDM_Lambda",            &RS.ELBDM_Lambda,            SID, TID, NonFatal, &RT.ELBDM_Lambda,             1, NonFatal );
#  endif
   LoadField( "ELBDM_Taylor3_Coeff",     &RS.ELBDM_Taylor3_Coeff,     SID, TID, NonFatal, &RT.ELBDM_Taylor3_Coeff,      1, NonFatal );
   LoadField( "ELBDM_Taylor3_Auto",      &RS.ELBDM_Taylor3_Auto,      SID, TID, NonFatal, &RT.ELBDM_Taylor3_Auto,       1, NonFatal );
#  endif

// fluid solvers in both HYDRO/ELBDM
   LoadField( "Flu_GPU_NPGroup",         &RS.Flu_GPU_NPGroup,         SID, TID, NonFatal, &RT.Flu_GPU_NPGroup,          1, NonFatal );
   LoadField( "GPU_NStream",             &RS.GPU_NStream,             SID, TID, NonFatal, &RT.GPU_NStream,              1, NonFatal );
   LoadField( "Opt__FixUp_Flux",         &RS.Opt__FixUp_Flux,         SID, TID, NonFatal, &RT.Opt__FixUp_Flux,          1, NonFatal );
#  ifdef MHD
   LoadField( "Opt__FixUp_Electric",     &RS.Opt__FixUp_Electric,     SID, TID, NonFatal, &RT.Opt__FixUp_Electric,      1, NonFatal );
#  endif
   LoadField( "Opt__FixUp_Restrict",     &RS.Opt__FixUp_Restrict,     SID, TID, NonFatal, &RT.Opt__FixUp_Restrict,      1, NonFatal );
   LoadField( "Opt__CorrAfterAllSync",   &RS.Opt__CorrAfterAllSync,   SID, TID, NonFatal, &RT.Opt__CorrAfterAllSync,    1, NonFatal );
   LoadField( "Opt__NormalizePassive",   &RS.Opt__NormalizePassive,   SID, TID, NonFatal, &RT.Opt__NormalizePassive,    1, NonFatal );
   LoadField( "NormalizePassive_NVar",   &RS.NormalizePassive_NVar,   SID, TID, NonFatal, &RT.NormalizePassive_NVar,    1, NonFatal );
   LoadField( "NormalizePassive_VarIdx",  RS.NormalizePassive_VarIdx, SID, TID, NonFatal,  RT.NormalizePassive_VarIdx, NP, NonFatal );
   LoadField( "Opt__OverlapMPI",         &RS.Opt__OverlapMPI,         SID, TID, NonFatal, &RT.Opt__OverlapMPI,          1, NonFatal );
   LoadField( "Opt__ResetFluid",         &RS.Opt__ResetFluid,         SID, TID, NonFatal, &RT.Opt__ResetFluid,          1, NonFatal );
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   LoadField( "MinDens",                 &RS.MinDens,                 SID, TID, NonFatal, &RT.MinDens,                  1, NonFatal );
#  endif
#  if ( MODEL == HYDRO )
   LoadField( "MinPres",                 &RS.MinPres,                 SID, TID, NonFatal, &RT.MinPres,                  1, NonFatal );
   LoadField( "MinEint",                 &RS.MinEint,                 SID, TID, NonFatal, &RT.MinEint,                  1, NonFatal );
   LoadField( "MinTemp",                 &RS.MinTemp,                 SID, TID, NonFatal, &RT.MinTemp,                  1, NonFatal );
   LoadField( "Opt__LastResortFloor",    &RS.Opt__LastResortFloor,    SID, TID, NonFatal, &RT.Opt__LastResortFloor,     1, NonFatal );
   LoadField( "JeansMinPres",            &RS.JeansMinPres,            SID, TID, NonFatal, &RT.JeansMinPres,             1, NonFatal );
   LoadField( "JeansMinPres_Level",      &RS.JeansMinPres_Level,      SID, TID, NonFatal, &RT.JeansMinPres_Level,       1, NonFatal );
   LoadField( "JeansMinPres_NCell",      &RS.JeansMinPres_NCell,      SID, TID, NonFatal, &RT.JeansMinPres_NCell,       1, NonFatal );
#  endif
#  ifdef DUAL_ENERGY
   LoadField( "DualEnergySwitch",        &RS.DualEnergySwitch,        SID, TID, NonFatal, &RT.DualEnergySwitch,         1, NonFatal );
#  endif

// self-gravity
#  ifdef GRAVITY
   LoadField( "NewtonG",                 &RS.NewtonG,                 SID, TID, NonFatal, &RT.NewtonG,                  1, NonFatal );
#  if   ( POT_SCHEME == SOR )
   LoadField( "SOR_Omega",               &RS.SOR_Omega,               SID, TID, NonFatal, &RT.SOR_Omega,                1, NonFatal );
   LoadField( "SOR_MaxIter",             &RS.SOR_MaxIter,             SID, TID, NonFatal, &RT.SOR_MaxIter,              1, NonFatal );
   LoadField( "SOR_MinIter",             &RS.SOR_MinIter,             SID, TID, NonFatal, &RT.SOR_MinIter,              1, NonFatal );
#  elif ( POT_SCHEME == MG )
   LoadField( "MG_MaxIter",              &RS.MG_MaxIter,              SID, TID, NonFatal, &RT.MG_MaxIter,               1, NonFatal );
   LoadField( "MG_NPreSmooth",           &RS.MG_NPreSmooth,           SID, TID, NonFatal, &RT.MG_NPreSmooth,            1, NonFatal );
   LoadField( "MG_NPostSmooth",          &RS.MG_NPostSmooth,          SID, TID, NonFatal, &RT.MG_NPostSmooth,           1, NonFatal );
   LoadField( "MG_ToleratedError",       &RS.MG_ToleratedError,       SID, TID, NonFatal, &RT.MG_ToleratedError,        1, NonFatal );
#  endif
   LoadField( "Pot_GPU_NPGroup",         &RS.Pot_GPU_NPGroup,         SID, TID, NonFatal, &RT.Pot_GPU_NPGroup,          1, NonFatal );
   LoadField( "Opt__GraP5Gradient",      &RS.Opt__GraP5Gradient,      SID, TID, NonFatal, &RT.Opt__GraP5Gradient,       1, NonFatal );
   LoadField( "Opt__SelfGravity",        &RS.Opt__SelfGravity,        SID, TID, NonFatal, &RT.Opt__SelfGravity,         1, NonFatal );
   LoadField( "Opt__ExtAcc",             &RS.Opt__ExtAcc,             SID, TID, NonFatal, &RT.Opt__ExtAcc,              1, NonFatal );
   LoadField( "Opt__ExtPot",             &RS.Opt__ExtPot,             SID, TID, NonFatal, &RT.Opt__ExtPot,              1, NonFatal );
   LoadField( "ExtPotTable_Name",        &RS.ExtPotTable_Name,        SID, TID, NonFatal,  RT.ExtPotTable_Name,         1, NonFatal );
   LoadField( "ExtPotTable_NPoint",       RS.ExtPotTable_NPoint,      SID, TID, NonFatal,  RT.ExtPotTable_NPoint,       3, NonFatal );
   LoadField( "ExtPotTable_dh",          &RS.ExtPotTable_dh,          SID, TID, NonFatal, &RT.ExtPotTable_dh,           1, NonFatal );
   LoadField( "ExtPotTable_EdgeL",        RS.ExtPotTable_EdgeL,       SID, TID, NonFatal,  RT.ExtPotTable_EdgeL,        3, NonFatal );
   LoadField( "ExtPotTable_Float8",      &RS.ExtPotTable_Float8,      SID, TID, NonFatal, &RT.ExtPotTable_Float8,       1, NonFatal );
   LoadField( "Opt__GravityExtraMass",   &RS.Opt__GravityExtraMass,   SID, TID, NonFatal, &RT.Opt__GravityExtraMass,    1, NonFatal );
#  endif

// source terms
   LoadField( "Src_Deleptonization",     &RS.Src_Deleptonization,     SID, TID, NonFatal, &RT.Src_Deleptonization,      1, NonFatal );
   LoadField( "Src_User",                &RS.Src_User,                SID, TID, NonFatal, &RT.Src_User,                 1, NonFatal );
   LoadField( "Src_GPU_NPGroup",         &RS.Src_GPU_NPGroup,         SID, TID, NonFatal, &RT.Src_GPU_NPGroup,          1, NonFatal );

// Grackle
#  ifdef SUPPORT_GRACKLE
   LoadField( "Grackle_Activate",        &RS.Grackle_Activate,        SID, TID, NonFatal, &RT.Grackle_Activate,         1, NonFatal );
   LoadField( "Grackle_Verbose",         &RS.Grackle_Verbose,         SID, TID, NonFatal, &RT.Grackle_Verbose,          1, NonFatal );
   LoadField( "Grackle_Cooling",         &RS.Grackle_Cooling,         SID, TID, NonFatal, &RT.Grackle_Cooling,          1, NonFatal );
   LoadField( "Grackle_Primordial",      &RS.Grackle_Primordial,      SID, TID, NonFatal, &RT.Grackle_Primordial,       1, NonFatal );
   LoadField( "Grackle_Metal",           &RS.Grackle_Metal,           SID, TID, NonFatal, &RT.Grackle_Metal,            1, NonFatal );
   LoadField( "Grackle_UV",              &RS.Grackle_UV,              SID, TID, NonFatal, &RT.Grackle_UV,               1, NonFatal );
   LoadField( "Grackle_CMB_Floor",       &RS.Grackle_CMB_Floor,       SID, TID, NonFatal, &RT.Grackle_CMB_Floor,        1, NonFatal );
   LoadField( "Grackle_PE_Heating",      &RS.Grackle_PE_Heating,      SID, TID, NonFatal, &RT.Grackle_PE_Heating,       1, NonFatal );
   LoadField( "Grackle_PE_HeatingRate",  &RS.Grackle_PE_HeatingRate,  SID, TID, NonFatal, &RT.Grackle_PE_HeatingRate,   1, NonFatal );
   LoadField( "Grackle_CloudyTable",     &RS.Grackle_CloudyTable,     SID, TID, NonFatal,  RT.Grackle_CloudyTable,      1, NonFatal );
   LoadField( "Grackle_ThreeBodyRate",   &RS.Grackle_ThreeBodyRate,   SID, TID, NonFatal, &RT.Grackle_ThreeBodyRate,    1, NonFatal );
   LoadField( "Grackle_CIE_Cooling",     &RS.Grackle_CIE_Cooling,     SID, TID, NonFatal, &RT.Grackle_CIE_Cooling,      1, NonFatal );
   LoadField( "Grackle_H2_OpaApprox",    &RS.Grackle_H2_OpaApprox,    SID, TID, NonFatal, &RT.Grackle_H2_OpaApprox,     1, NonFatal );
   LoadField( "Che_GPU_NPGroup",         &RS.Che_GPU_NPGroup,         SID, TID, NonFatal, &RT.Che_GPU_NPGroup,          1, NonFatal );
#  endif

// star formation
#  ifdef STAR_FORMATION
   LoadField( "SF_CreateStar_Scheme",       &RS.SF_CreateStar_Scheme,       SID, TID, NonFatal, &RT.SF_CreateStar_Scheme,       1, NonFatal );
   LoadField( "SF_CreateStar_RSeed",        &RS.SF_CreateStar_RSeed,        SID, TID, NonFatal, &RT.SF_CreateStar_RSeed,        1, NonFatal );
   LoadField( "SF_CreateStar_DetRandom",    &RS.SF_CreateStar_DetRandom,    SID, TID, NonFatal, &RT.SF_CreateStar_DetRandom,    1, NonFatal );
   LoadField( "SF_CreateStar_MinLevel",     &RS.SF_CreateStar_MinLevel,     SID, TID, NonFatal, &RT.SF_CreateStar_MinLevel,     1, NonFatal );
   LoadField( "SF_CreateStar_MinGasDens",   &RS.SF_CreateStar_MinGasDens,   SID, TID, NonFatal, &RT.SF_CreateStar_MinGasDens,   1, NonFatal );
   LoadField( "SF_CreateStar_MassEff",      &RS.SF_CreateStar_MassEff,      SID, TID, NonFatal, &RT.SF_CreateStar_MassEff,      1, NonFatal );
   LoadField( "SF_CreateStar_MinStarMass",  &RS.SF_CreateStar_MinStarMass,  SID, TID, NonFatal, &RT.SF_CreateStar_MinStarMass,  1, NonFatal );
   LoadField( "SF_CreateStar_MaxStarMFrac", &RS.SF_CreateStar_MaxStarMFrac, SID, TID, NonFatal, &RT.SF_CreateStar_MaxStarMFrac, 1, NonFatal );
#  endif

// initialization
   LoadField( "Opt__Init",               &RS.Opt__Init,               SID, TID, NonFatal, &RT.Opt__Init,                1, NonFatal );
   LoadField( "RestartLoadNRank",        &RS.RestartLoadNRank,        SID, TID, NonFatal, &RT.RestartLoadNRank,         1, NonFatal );
   LoadField( "Opt__RestartReset",       &RS.Opt__RestartReset,       SID, TID, NonFatal, &RT.Opt__RestartReset,        1, NonFatal );
   LoadField( "Opt__UM_IC_Level",        &RS.Opt__UM_IC_Level,        SID, TID, NonFatal, &RT.Opt__UM_IC_Level,         1, NonFatal );
   LoadField( "Opt__UM_IC_NVar",         &RS.Opt__UM_IC_NVar,         SID, TID, NonFatal, &RT.Opt__UM_IC_NVar,          1, NonFatal );
   LoadField( "Opt__UM_IC_Format",       &RS.Opt__UM_IC_Format,       SID, TID, NonFatal, &RT.Opt__UM_IC_Format,        1, NonFatal );
   LoadField( "Opt__UM_IC_Downgrade",    &RS.Opt__UM_IC_Downgrade,    SID, TID, NonFatal, &RT.Opt__UM_IC_Downgrade,     1, NonFatal );
   LoadField( "Opt__UM_IC_Refine",       &RS.Opt__UM_IC_Refine,       SID, TID, NonFatal, &RT.Opt__UM_IC_Refine,        1, NonFatal );
   LoadField( "Opt__UM_IC_LoadNRank",    &RS.Opt__UM_IC_LoadNRank,    SID, TID, NonFatal, &RT.Opt__UM_IC_LoadNRank,     1, NonFatal );
   LoadField( "Opt__InitRestrict",       &RS.Opt__InitRestrict,       SID, TID, NonFatal, &RT.Opt__InitRestrict,        1, NonFatal );
   LoadField( "Opt__InitGridWithOMP",    &RS.Opt__InitGridWithOMP,    SID, TID, NonFatal, &RT.Opt__InitGridWithOMP,     1, NonFatal );
   LoadField( "Opt__GPUID_Select",       &RS.Opt__GPUID_Select,       SID, TID, NonFatal, &RT.Opt__GPUID_Select,        1, NonFatal );
   LoadField( "Init_Subsampling_NCell",  &RS.Init_Subsampling_NCell,  SID, TID, NonFatal, &RT.Init_Subsampling_NCell,   1, NonFatal );

// interpolation schemes
   LoadField( "Opt__Int_Time",           &RS.Opt__Int_Time,           SID, TID, NonFatal, &RT.Opt__Int_Time,            1, NonFatal );
#  if ( MODEL == ELBDM )
   LoadField( "Opt__Int_Phase",          &RS.Opt__Int_Phase,          SID, TID, NonFatal, &RT.Opt__Int_Phase,           1, NonFatal );
#  endif
   LoadField( "Opt__Flu_IntScheme",      &RS.Opt__Flu_IntScheme,      SID, TID, NonFatal, &RT.Opt__Flu_IntScheme,       1, NonFatal );
   LoadField( "Opt__RefFlu_IntScheme",   &RS.Opt__RefFlu_IntScheme,   SID, TID, NonFatal, &RT.Opt__RefFlu_IntScheme,    1, NonFatal );
#  ifdef MHD
   LoadField( "Opt__Mag_IntScheme",      &RS.Opt__Mag_IntScheme,      SID, TID, NonFatal, &RT.Opt__Mag_IntScheme,       1, NonFatal );
   LoadField( "Opt__RefMag_IntScheme",   &RS.Opt__RefMag_IntScheme,   SID, TID, NonFatal, &RT.Opt__RefMag_IntScheme,    1, NonFatal );
#  endif
#  ifdef GRAVITY
   LoadField( "Opt__Pot_IntScheme",      &RS.Opt__Pot_IntScheme,      SID, TID, NonFatal, &RT.Opt__Pot_IntScheme,       1, NonFatal );
   LoadField( "Opt__Rho_IntScheme",      &RS.Opt__Rho_IntScheme,      SID, TID, NonFatal, &RT.Opt__Rho_IntScheme,       1, NonFatal );
   LoadField( "Opt__Gra_IntScheme",      &RS.Opt__Gra_IntScheme,      SID, TID, NonFatal, &RT.Opt__Gra_IntScheme,       1, NonFatal );
   LoadField( "Opt__RefPot_IntScheme",   &RS.Opt__RefPot_IntScheme,   SID, TID, NonFatal, &RT.Opt__RefPot_IntScheme,    1, NonFatal );
#  endif
   LoadField( "IntMonoCoeff",            &RS.IntMonoCoeff,            SID, TID, NonFatal, &RT.IntMonoCoeff,             1, NonFatal );
   LoadField( "IntOppSign0thOrder",      &RS.IntOppSign0thOrder,      SID, TID, NonFatal, &RT.IntOppSign0thOrder,       1, NonFatal );

// data dump
   LoadField( "Opt__Output_Total",       &RS.Opt__Output_Total,       SID, TID, NonFatal, &RT.Opt__Output_Total,        1, NonFatal );
   LoadField( "Opt__Output_Part",        &RS.Opt__Output_Part,        SID, TID, NonFatal, &RT.Opt__Output_Part,         1, NonFatal );
   LoadField( "Opt__Output_User",        &RS.Opt__Output_User,        SID, TID, NonFatal, &RT.Opt__Output_User,         1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Output_ParText",     &RS.Opt__Output_ParText,     SID, TID, NonFatal, &RT.Opt__Output_ParText,      1, NonFatal );
#  endif
   LoadField( "Opt__Output_BasePS",      &RS.Opt__Output_BasePS,      SID, TID, NonFatal, &RT.Opt__Output_BasePS,       1, NonFatal );
   if ( OPT__OUTPUT_PART )
   LoadField( "Opt__Output_Base",        &RS.Opt__Output_Base,        SID, TID, NonFatal, &RT.Opt__Output_Base,         1, NonFatal );
#  ifdef GRAVITY
   LoadField( "Opt__Output_Pot",         &RS.Opt__Output_Pot,         SID, TID, NonFatal, &RT.Opt__Output_Pot,          1, NonFatal );
#  endif
#  ifdef PARTICLE
   LoadField( "Opt__Output_ParDens",     &RS.Opt__Output_ParDens,     SID, TID, NonFatal, &RT.Opt__Output_ParDens,      1, NonFatal );
#  endif
#  ifdef MHD
   LoadField( "Opt__Output_CC_Mag",      &RS.Opt__Output_CC_Mag,      SID, TID, NonFatal, &RT.Opt__Output_CC_Mag,       1, NonFatal );
#  endif
#  if ( MODEL == HYDRO )
   LoadField( "Opt__Output_Pres",        &RS.Opt__Output_Pres,        SID, TID, NonFatal, &RT.Opt__Output_Pres,         1, NonFatal );
   LoadField( "Opt__Output_Temp",        &RS.Opt__Output_Temp,        SID, TID, NonFatal, &RT.Opt__Output_Temp,         1, NonFatal );
   LoadField( "Opt__Output_Cs",          &RS.Opt__Output_Cs,          SID, TID, NonFatal, &RT.Opt__Output_Cs,           1, NonFatal );
   LoadField( "Opt__Output_DivVel",      &RS.Opt__Output_DivVel,      SID, TID, NonFatal, &RT.Opt__Output_DivVel,       1, NonFatal );
   LoadField( "Opt__Output_Mach",        &RS.Opt__Output_Mach,        SID, TID, NonFatal, &RT.Opt__Output_Mach,         1, NonFatal );
#  ifdef MHD
   LoadField( "Opt__Output_DivMag",      &RS.Opt__Output_DivMag,      SID, TID, NonFatal, &RT.Opt__Output_DivMag,       1, NonFatal );
#  endif
#  endif // #if ( MODEL == HYDRO )
   LoadField( "Opt__Output_UserField",   &RS.Opt__Output_UserField,   SID, TID, NonFatal, &RT.Opt__Output_UserField,    1, NonFatal );
#  ifdef PARTICLE
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PAR_TEXT ) {
#  else
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS ) {
#  endif
   LoadField( "Opt__Output_Mode",        &RS.Opt__Output_Mode,        SID, TID, NonFatal, &RT.Opt__Output_Mode,         1, NonFatal );
   LoadField( "Opt__Output_Step",        &RS.Opt__Output_Step,        SID, TID, NonFatal, &RT.Opt__Output_Step,         1, NonFatal );
   LoadField( "Opt__Output_Dt",          &RS.Opt__Output_Dt,          SID, TID, NonFatal, &RT.Opt__Output_Dt,           1, NonFatal );
   }
   if ( OPT__OUTPUT_PART ) {
   LoadField( "Output_PartX",            &RS.Output_PartX,            SID, TID, NonFatal, &RT.Output_PartX,             1, NonFatal );
   LoadField( "Output_PartY",            &RS.Output_PartY,            SID, TID, NonFatal, &RT.Output_PartY,             1, NonFatal );
   LoadField( "Output_PartZ",            &RS.Output_PartZ,            SID, TID, NonFatal, &RT.Output_PartZ,             1, NonFatal );
   }
   LoadField( "InitDumpID",              &RS.InitDumpID,              SID, TID, NonFatal, &RT.InitDumpID,               1, NonFatal );

// miscellaneous
   LoadField( "Opt__Verbose",            &RS.Opt__Verbose,            SID, TID, NonFatal, &RT.Opt__Verbose,             1, NonFatal );
   LoadField( "Opt__TimingBarrier",      &RS.Opt__TimingBarrier,      SID, TID, NonFatal, &RT.Opt__TimingBarrier,       1, NonFatal );
   LoadField( "Opt__TimingBalance",      &RS.Opt__TimingBalance,      SID, TID, NonFatal, &RT.Opt__TimingBalance,       1, NonFatal );
   LoadField( "Opt__TimingMPI",          &RS.Opt__TimingMPI,          SID, TID, NonFatal, &RT.Opt__TimingMPI,           1, NonFatal );
   LoadField( "Opt__RecordNote",         &RS.Opt__RecordNote,         SID, TID, NonFatal, &RT.Opt__RecordNote,          1, NonFatal );
   LoadField( "Opt__RecordUnphy",        &RS.Opt__RecordUnphy,        SID, TID, NonFatal, &RT.Opt__RecordUnphy,         1, NonFatal );
   LoadField( "Opt__RecordMemory",       &RS.Opt__RecordMemory,       SID, TID, NonFatal, &RT.Opt__RecordMemory,        1, NonFatal );
   LoadField( "Opt__RecordPerformance",  &RS.Opt__RecordPerformance,  SID, TID, NonFatal, &RT.Opt__RecordPerformance,   1, NonFatal );
   LoadField( "Opt__ManualControl",      &RS.Opt__ManualControl,      SID, TID, NonFatal, &RT.Opt__ManualControl,       1, NonFatal );
   LoadField( "Opt__RecordUser",         &RS.Opt__RecordUser,         SID, TID, NonFatal, &RT.Opt__RecordUser,          1, NonFatal );
   LoadField( "Opt__OptimizeAggressive", &RS.Opt__OptimizeAggressive, SID, TID, NonFatal, &RT.Opt__OptimizeAggressive,  1, NonFatal );

// simulation checks
   LoadField( "Opt__Ck_Refine",          &RS.Opt__Ck_Refine,          SID, TID, NonFatal, &RT.Opt__Ck_Refine,           1, NonFatal );
   LoadField( "Opt__Ck_ProperNesting",   &RS.Opt__Ck_ProperNesting,   SID, TID, NonFatal, &RT.Opt__Ck_ProperNesting,    1, NonFatal );
   LoadField( "Opt__Ck_Conservation",    &RS.Opt__Ck_Conservation,    SID, TID, NonFatal, &RT.Opt__Ck_Conservation,     1, NonFatal );
   LoadField( "Opt__Ck_NormPassive",     &RS.Opt__Ck_NormPassive,     SID, TID, NonFatal, &RT.Opt__Ck_NormPassive,      1, NonFatal );
   LoadField( "Opt__Ck_Restrict",        &RS.Opt__Ck_Restrict,        SID, TID, NonFatal, &RT.Opt__Ck_Restrict,         1, NonFatal );
   LoadField( "Opt__Ck_Finite",          &RS.Opt__Ck_Finite,          SID, TID, NonFatal, &RT.Opt__Ck_Finite,           1, NonFatal );
   LoadField( "Opt__Ck_PatchAllocate",   &RS.Opt__Ck_PatchAllocate,   SID, TID, NonFatal, &RT.Opt__Ck_PatchAllocate,    1, NonFatal );
   LoadField( "Opt__Ck_FluxAllocate",    &RS.Opt__Ck_FluxAllocate,    SID, TID, NonFatal, &RT.Opt__Ck_FluxAllocate,     1, NonFatal );
#  if ( MODEL == HYDRO )
   LoadField( "Opt__Ck_Negative",        &RS.Opt__Ck_Negative,        SID, TID, NonFatal, &RT.Opt__Ck_Negative,         1, NonFatal );
#  endif
   LoadField( "Opt__Ck_MemFree",         &RS.Opt__Ck_MemFree,         SID, TID, NonFatal, &RT.Opt__Ck_MemFree,          1, NonFatal );
#  ifdef PARTICLE
   LoadField( "Opt__Ck_Particle",        &RS.Opt__Ck_Particle,        SID, TID, NonFatal, &RT.Opt__Ck_Particle,         1, NonFatal );
#  endif
#  ifdef MHD
   LoadField( "Opt__Ck_InterfaceB",      &RS.Opt__Ck_InterfaceB,      SID, TID, NonFatal, &RT.Opt__Ck_InterfaceB,       1, NonFatal );
   LoadField( "Opt__Ck_DivergenceB",     &RS.Opt__Ck_DivergenceB,     SID, TID, NonFatal, &RT.Opt__Ck_DivergenceB,      1, NonFatal );
#  endif


// flag tables
#  if   ( MODEL == HYDRO )
   const bool Opt__FlagLohner = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP );
#  elif ( MODEL == ELBDM )
   const bool Opt__FlagLohner = OPT__FLAG_LOHNER_DENS;
#  else
#  error : ERROR : unsupported MODEL !!
#  endif

// initialize as -1 (to work with NLvRestart < NLEVEL)
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      RS.FlagTable_Rho         [lv]    = -1.0;
      RS.FlagTable_RhoGradient [lv]    = -1.0;

      for (int t=0; t<4; t++)
      RS.FlagTable_Lohner      [lv][t] = -1.0;

      RS.FlagTable_User        [lv].p   = malloc( OPT__FLAG_USER_NUM*sizeof(double) );
      RS.FlagTable_User        [lv].len = OPT__FLAG_USER_NUM;
      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      ( (double *) RS.FlagTable_User[lv].p )[t] = -1.0;

#     if   ( MODEL == HYDRO )
      RS.FlagTable_PresGradient[lv]    = -1.0;
      RS.FlagTable_Vorticity   [lv]    = -1.0;
      RS.FlagTable_Jeans       [lv]    = -1.0;
#     ifdef MHD
      RS.FlagTable_Current     [lv]    = -1.0;
#     endif

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++)
      RS.FlagTable_EngyDensity [lv][t] = -1.0;
#     endif

#     ifdef PARTICLE
      RS.FlagTable_NParPatch   [lv]    = -1;
      RS.FlagTable_NParCell    [lv]    = -1;
      RS.FlagTable_ParMassCell [lv]    = -1.0;
#     endif
   }

   if ( OPT__FLAG_RHO )
   LoadField( "FlagTable_Rho",            RS.FlagTable_Rho,           SID, TID, NonFatal,  RT.FlagTable_Rho,           N1, NonFatal );

   if ( OPT__FLAG_RHO_GRADIENT )
   LoadField( "FlagTable_RhoGradient",    RS.FlagTable_RhoGradient,   SID, TID, NonFatal,  RT.FlagTable_RhoGradient,   N1, NonFatal );

   if ( Opt__FlagLohner ) {
   LoadField( "FlagTable_Lohner",         RS.FlagTable_Lohner,        SID, TID, NonFatal,  NullPtr,                    -1, NonFatal );

   for (int lv=0; lv<MAX_LEVEL; lv++)
   for (int t=0; t<4; t++)
   {
      if ( RS.FlagTable_Lohner[lv][t] != RT.FlagTable_Lohner[lv][t] )
         Aux_Message( stderr, "WARNING : \"%s[%d][%d]\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                       "FlagTable_Lohner", lv, t, RS.FlagTable_Lohner[lv][t],  RT.FlagTable_Lohner[lv][t] );
   }}

   if ( OPT__FLAG_USER ) {
   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
      char Key[MAX_STRING];
      sprintf( Key, "FlagTable_User_Lv%02d", lv );

      LoadField( Key,                    &RS.FlagTable_User[lv],      SID, TID, NonFatal,  NullPtr,                    -1, NonFatal );

      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      if (  ( (double *) RS.FlagTable_User[lv].p )[t] != ( (double *) RT.FlagTable_User[lv].p )[t]  )
         Aux_Message( stderr, "WARNING : \"%s[%d][%d]\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                      "FlagTable_User", lv, t, ( (double *) RS.FlagTable_User[lv].p )[t],  ( (double *) RT.FlagTable_User[lv].p )[t] );
   }}

#  if   ( MODEL == HYDRO )
   if ( OPT__FLAG_PRES_GRADIENT )
   LoadField( "FlagTable_PresGradient",   RS.FlagTable_PresGradient,  SID, TID, NonFatal,  RT.FlagTable_PresGradient,  N1, NonFatal );

   if ( OPT__FLAG_VORTICITY )
   LoadField( "FlagTable_Vorticity",      RS.FlagTable_Vorticity,     SID, TID, NonFatal,  RT.FlagTable_Vorticity,     N1, NonFatal );

   if ( OPT__FLAG_JEANS )
   LoadField( "FlagTable_Jeans",          RS.FlagTable_Jeans,         SID, TID, NonFatal,  RT.FlagTable_Jeans,         N1, NonFatal );

#  ifdef MHD
   if ( OPT__FLAG_CURRENT )
   LoadField( "FlagTable_Current",        RS.FlagTable_Current,       SID, TID, NonFatal,  RT.FlagTable_Current,       N1, NonFatal );
#  endif

#  elif ( MODEL == ELBDM )
   if ( OPT__FLAG_ENGY_DENSITY ) {
   LoadField( "FlagTable_EngyDensity",    RS.FlagTable_EngyDensity,   SID, TID, NonFatal,  NullPtr,                    -1, NonFatal );

   for (int lv=0; lv<MAX_LEVEL; lv++)
   for (int t=0; t<2; t++)
   {
      if ( RS.FlagTable_EngyDensity[lv][t] != RT.FlagTable_EngyDensity[lv][t] )
         Aux_Message( stderr, "WARNING : \"%s[%d][%d]\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                       "FlagTable_EngyDensity", lv, t, RS.FlagTable_EngyDensity[lv][t],  RT.FlagTable_EngyDensity[lv][t] );
   }}
#  endif

#  ifdef PARTICLE
   if ( OPT__FLAG_NPAR_PATCH )
   LoadField( "FlagTable_NParPatch",      RS.FlagTable_NParPatch,     SID, TID, NonFatal,  RT.FlagTable_NParPatch,     N1, NonFatal );

   if ( OPT__FLAG_NPAR_CELL )
   LoadField( "FlagTable_NParCell",       RS.FlagTable_NParCell,      SID, TID, NonFatal,  RT.FlagTable_NParCell,      N1, NonFatal );

   if ( OPT__FLAG_PAR_MASS_CELL )
   LoadField( "FlagTable_ParMassCell",    RS.FlagTable_ParMassCell,   SID, TID, NonFatal,  RT.FlagTable_ParMassCell,   N1, NonFatal );
#  endif


// 5. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

// free memory
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      free( RS.FlagTable_User[lv].p );
      free( RT.FlagTable_User[lv].p );
   }

} // FUNCTION : Check_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  ResetParameter
// Description :  Reset simulation parameters from the restart file
//
// Note        :  1. Replace some runtime parameters by the values loaded from the restart file
//                2. This function must be called by ALL ranks
//
// Parameter   :  FileName : Restart file name
//                EndT     : Runtime parameter "END_T"
//                EndStep  : Runtime parameter "END_STEP"
//-------------------------------------------------------------------------------------------------------
void ResetParameter( const char *FileName, double *EndT, long *EndStep )
{

   const bool    Fatal = true;
   const bool NonFatal = false;
   const int *NullPtr  = NULL;

   herr_t Status;


// 1. open the HDF5 file
   const hid_t FID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );  // file ID

   if ( FID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 2. load the InputPar dataset and datatype
   const hid_t SID = H5Dopen( FID, "Info/InputPara", H5P_DEFAULT ); // dataset ID

   if ( SID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );

   const hid_t TID = H5Dget_type( SID );                                // datatype ID


// 3. reset the runtime parameters if required
   InputPara_t RS;   // RS = ReStart

   if ( *EndT < 0.0 )
   {
      LoadField( "EndT",      EndT,    SID, TID, Fatal, NullPtr, -1, NonFatal );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %14.7e\n", "END_T", *EndT );
   }

   if ( *EndStep < 0 )
   {
      LoadField( "EndStep",   EndStep, SID, TID, Fatal, NullPtr, -1, NonFatal );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %ld\n", "END_STEP", *EndStep );
   }

// 4. close all objects
   Status = H5Tclose( TID );
   Status = H5Dclose( SID );
   Status = H5Fclose( FID );

} // FUNCTION : ResetParameter



#endif // #ifdef SUPPORT_HDF5
