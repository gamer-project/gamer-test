#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif
#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

void Init_ByRestart_v1( const char FileName[] );
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, int &NLv_Restart,
                                bool &LoadPot, bool &LoadParDens, bool &LoadCCMag,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter );
void CompareVar( const char *VarName, const bool   RestartVar, const bool   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const int    RestartVar, const int    RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const long   RestartVar, const long   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const real   RestartVar, const real   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByRestart
// Description :  Reload a previous output as the initial condition
//
// Note        :  1. This function will always load the file named "RESTART"
//                   --> You can just make a symbolic link named RESTART to the file you want to use as the
//                       initial condition
//
//                2. This function will invoke "Init_ByRestart_HDF5" automatically if the restart file
//                   is in the HDF5 format
//
//                3. This function will invoke "Init_ByRestart_v1" automatically if the restart file
//                   is in a simple binary format in version 1 (i.e., FormatVersion < 2000)
//-------------------------------------------------------------------------------------------------------
void Init_ByRestart()
{

   const char FileName[] = "RESTART";


// load the HDF5 data
#  ifdef SUPPORT_HDF5
   if (  Aux_CheckFileExist(FileName)  &&  H5Fis_hdf5(FileName)  )
   {
      Init_ByRestart_HDF5( FileName );
      return;
   }
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if the restart file exists
   if ( !Aux_CheckFileExist(FileName)  &&  MPI_Rank == 0 )
      Aux_Error( ERROR_INFO, "restart file \"%s\" does not exist !!\n", FileName );

   MPI_Barrier( MPI_COMM_WORLD );


   FILE *File = fopen( FileName, "rb" );

// record the size of different data types
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


// a. load the information of data format
// =================================================================================================
   long FormatVersion, CheckCode;
   long HeaderSize_Format, HeaderSize_Makefile, HeaderSize_Constant, HeaderSize_Parameter, HeaderSize_SimuInfo;
   long HeaderOffset_Format, HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, HeaderOffset_SimuInfo;
   long HeaderSize_Total;


// load the file format
   fread( &FormatVersion, sizeof(long), 1, File );


// switch to version 1
   if ( FormatVersion < 2000 )
   {
      rewind( File );
      fclose( File );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Switch to version 1 ...\n" );

      Init_ByRestart_v1( FileName );
      return;
   }


// load the size of each header (in bytes)
   fread( &HeaderSize_Format,    sizeof(long), 1, File );
   fread( &HeaderSize_Makefile,  sizeof(long), 1, File );
   fread( &HeaderSize_Constant,  sizeof(long), 1, File );
   fread( &HeaderSize_Parameter, sizeof(long), 1, File );
   fread( &HeaderSize_SimuInfo,  sizeof(long), 1, File );
   fread( &CheckCode,            sizeof(long), 1, File );


// determine the offset of each header (in bytes)
   HeaderOffset_Format    = 0;    // it must be zero
   HeaderOffset_Makefile  = HeaderOffset_Format    + HeaderSize_Format;
   HeaderOffset_Constant  = HeaderOffset_Makefile  + HeaderSize_Makefile;
   HeaderOffset_Parameter = HeaderOffset_Constant  + HeaderSize_Constant;
   HeaderOffset_SimuInfo  = HeaderOffset_Parameter + HeaderSize_Parameter;

   HeaderSize_Total       = HeaderOffset_SimuInfo  + HeaderSize_SimuInfo;


// verify the input data format version
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "   The format version of the RESTART file = %ld\n", FormatVersion );

      if ( FormatVersion < 2000 )
         Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2000) !!\n" );

#     ifdef PARTICLE
      if ( FormatVersion < 2100 )
         Aux_Error( ERROR_INFO, "unsupported data format version for particles (only support version >= 2100) !!\n" );
#     endif

#     ifdef MHD
      if ( FormatVersion < 2210 )
         Aux_Error( ERROR_INFO, "unsupported data format version for MHD (only support version >= 2210) !!\n" );
#     endif
   }
   MPI_Barrier( MPI_COMM_WORLD );


// check if the size of different data types are consistent
   int size_bool_restart, size_int_restart, size_long_restart, size_real_restart, size_double_restart;

   fread( &size_bool_restart,   sizeof(int), 1, File );
   fread( &size_int_restart,    sizeof(int), 1, File );
   fread( &size_long_restart,   sizeof(int), 1, File );
   fread( &size_real_restart,   sizeof(int), 1, File );
   fread( &size_double_restart, sizeof(int), 1, File );

   if ( size_bool_restart != size_bool )
      Aux_Error( ERROR_INFO, "sizeof(bool) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_bool_restart, size_bool );

   if ( size_int_restart != size_int )
      Aux_Error( ERROR_INFO, "sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_int_restart, size_int );

   if ( size_long_restart != size_long )
      Aux_Error( ERROR_INFO, "sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_long_restart, size_long );

   if ( size_real_restart != size_real )
      Aux_Error( ERROR_INFO, "sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_real_restart, size_real );

   if ( size_double_restart != size_double )
      Aux_Error( ERROR_INFO, "sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                 size_double_restart, size_double );


// b. load all simulation parameters
// =================================================================================================
   int  NLv_Restart = NLEVEL;
   bool LoadPot     = false;
   bool LoadParDens = false;
   bool LoadCCMag   = false;

   Load_Parameter_After_2000( File, FormatVersion, NLv_Restart, LoadPot, LoadParDens, LoadCCMag,
                              HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter );


// set the rescale factor for different NLEVEL
   const int rescale = 1 << ( NLEVEL - NLv_Restart );
   if ( MPI_Rank == 0  &&  rescale != 1 )
      Aux_Message( stderr, "WARNING : rescale factor is set to %d\n", rescale );



// c. load the simulation information
// =================================================================================================
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ...\n" );

// verify the check code
   long checkcode;

   fseek( File, HeaderOffset_SimuInfo, SEEK_SET );

   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
      Aux_Error( ERROR_INFO, "incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
                 checkcode, CheckCode );


// load information necessary for restart
   int  NDataPatch_Total[NLv_Restart];
#  ifdef PARTICLE
   long FileOffset_Particle;
#  endif

   fread( &DumpID,                        sizeof(int),              1, File );
   fread( Time,                           sizeof(double), NLv_Restart, File );
   fread( &Step,                          sizeof(long),             1, File );
   fread( NPatchTotal,                    sizeof(int),    NLv_Restart, File );
   fread( NDataPatch_Total,               sizeof(int),    NLv_Restart, File );
   fread( AdvanceCounter,                 sizeof(long),   NLv_Restart, File );

#  ifdef GRAVITY
   fread( &AveDensity_Init,               sizeof(double),           1, File );
#  else
   fseek( File, sizeof(double), SEEK_CUR );
#  endif

#  ifdef PARTICLE
   fread( &amr->Par->NPar_Active_AllRank, sizeof(long),             1, File );
   fread( &FileOffset_Particle,           sizeof(long),             1, File );
#  else
   fseek( File, 2*sizeof(long), SEEK_CUR );
#  endif

   if ( FormatVersion >= 2130 )
   fread( dTime_AllLv,                    sizeof(double), NLv_Restart, File );


// set parameters in levels that do not exist in the input file
// --> assuming dTime_AllLv[] has been initialized as 0.0 properly
   for (int lv=NLv_Restart; lv<NLEVEL; lv++)
   {
      Time          [lv] = Time[0];
      NPatchTotal   [lv] = 0;
      AdvanceCounter[lv] = 0;
   }


// reset parameters for OPT__RESTART_RESET
   if ( OPT__RESTART_RESET )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
#        ifdef COMOVING
         Time          [lv] = A_INIT;
#        else
         Time          [lv] = 0.0;
#        endif
         AdvanceCounter[lv] = 0;
         dTime_AllLv   [lv] = 0.0;
      }

      Step            = 0;
#     ifdef GRAVITY
      AveDensity_Init = -1.0;    // set to an arbitrary negative value
#     endif
   }


// set Flu(Pot)SgTime
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


// set the next dump ID
   if ( INIT_DUMPID < 0 )  DumpID = ( OPT__RESTART_RESET ) ? 0 : DumpID+1;
   else                    DumpID = INIT_DUMPID;


// verify the size of the RESTART file
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Verifying the size of the RESTART file ...\n" );

   long ExpectSize, InputSize, PatchDataSize, DataSize[NLv_Restart];
   int  NGridVar = NCOMP_TOTAL;  // number of grid variables

#  ifdef GRAVITY
   if ( LoadPot )       NGridVar ++;
#  endif
#  ifdef PARTICLE
   if ( LoadParDens )   NGridVar ++;
#  endif
#  ifdef MHD
   if ( LoadCCMag )     NGridVar += NCOMP_MAG;
#  endif

   PatchDataSize  = CUBE(PS1)*NGridVar*sizeof(real);
#  ifdef MHD
   PatchDataSize += PS1P1*SQR(PS1)*NCOMP_MAG*sizeof(real);
#  endif
   ExpectSize     = HeaderSize_Total;

   for (int lv=0; lv<NLv_Restart; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += (long)NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += (long)NDataPatch_Total[lv]*PatchDataSize;

      ExpectSize   += (long)DataSize[lv];
   }

#  ifdef PARTICLE
   for (int lv=0; lv<NLv_Restart; lv++)
   {
//    2 = NPar + starting particle index stored in each leaf patch
      const long ParInfoSize = (long)NDataPatch_Total[lv]*2*sizeof(long);

      DataSize[lv] += ParInfoSize;
      ExpectSize   += ParInfoSize;
   }

   ExpectSize += (long)PAR_NATT_STORED*amr->Par->NPar_Active_AllRank*sizeof(real);
#  endif

   fseek( File, 0, SEEK_END );
   InputSize = ftell( File );

   if ( InputSize != ExpectSize  &&  MPI_Rank == 0 )
      Aux_Error( ERROR_INFO, "size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
                 FileName, InputSize, ExpectSize );

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      Verifying the size of the RESTART file ... passed\n" );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation information ... done\n" );



// d. load the simulation grid data
// =================================================================================================
   int   Load_Cr_and_Son[4];
   int  *LoadCorner  = Load_Cr_and_Son;
   int  *LoadSon     = Load_Cr_and_Son + 3;
#  ifdef PARTICLE
   long  Load_NPar_and_GParID[2];
   long *Load_NPar   = Load_NPar_and_GParID;
   long *Load_GParID = Load_NPar_and_GParID + 1;

   int   MaxNParInOnePatch = 0;
   long  NParThisRank      = 0;
#  endif

// d0. set the load-balance cut points
#  ifdef LOAD_BALANCE
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting load-balance cut points ...\n" );

   const bool InputLBIdx0AndLoad_Yes = true;
   long   *LBIdx0_AllRank = NULL;
   double *Load_AllRank   = NULL;

   if ( MPI_Rank == 0 )
   {
      File = fopen( FileName, "rb" );
      fseek( File, HeaderSize_Total, SEEK_SET );
   }

   for (int lv=0; lv<NLv_Restart; lv++)
   {
//    d0-1. construct the LBIdx0_AllRank and Load_AllRank lists at rank 0
      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = new long   [ NPatchTotal[lv] / 8 ];
         Load_AllRank   = new double [ NPatchTotal[lv] / 8 ];

//       LBIdx0_AllRank
         for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
         {
//          load the corner and son of this patch
            fread( Load_Cr_and_Son, sizeof(int), 4, File );

//          only store the minimum LBIdx in each patch group
            if ( LoadPID%8 == 0 )
            {
               for (int d=0; d<3; d++)    LoadCorner[d] *= rescale;

               const int t = LoadPID / 8;

               LBIdx0_AllRank[t]  = LB_Corner2Index( lv, LoadCorner, CHECK_ON );
               LBIdx0_AllRank[t] -= LBIdx0_AllRank[t] % 8;
            }

            if ( *LoadSon == -1 )
            {
//             for particles, skip NPar and GParID as well
#              ifdef PARTICLE
               fseek( File, PatchDataSize+2*sizeof(long), SEEK_CUR );
#              else
               fseek( File, PatchDataSize,                SEEK_CUR );
#              endif
            }
         } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

//       Load_AllRank (assuming all patches have the same weighting == 1.0)
         for (int t=0; t<NPatchTotal[lv]/8; t++)   Load_AllRank[t] = 8.0;  // 8 patches per patch group
      } // if ( MPI_Rank == 0 )

//    d0-2. set the cut points
//    --> do NOT consider load-balance weighting of particles since at this point we don't have that information
      const double ParWeight_Zero = 0.0;
      LB_SetCutPoint( lv, NPatchTotal[lv]/8, amr->LB->CutPoint[lv], InputLBIdx0AndLoad_Yes, LBIdx0_AllRank,
                      Load_AllRank, ParWeight_Zero );

      if ( MPI_Rank == 0 )
      {
         delete [] LBIdx0_AllRank;
         delete [] Load_AllRank;
      }
   } // for (int lv=0; lv<NLv_Restart; lv++)

   if ( MPI_Rank == 0 )    fclose( File );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting load-balance cut points ... done\n" );
#  endif // #ifdef LOAD_BALANCE


// begin to load data
   long Offset = HeaderSize_Total;
   int  PID;
#  ifndef LOAD_BALANCE
   int TargetRange_Min[3], TargetRange_Max[3];
#  endif

   for (int lv=0; lv<NLv_Restart; lv++)
   {
      for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)
      {
         if ( MPI_Rank == 0 )
            Aux_Message( stdout, "   Loading grid data at level %2d, MPI ranks %4d -- %4d ... ",
                         lv, TRanks, MIN(TRanks+RESTART_LOAD_NRANK-1, MPI_NRank-1) );

         if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )
         {
//          d1. set the range of the target sub-domain
#           ifndef LOAD_BALANCE
            for (int d=0; d<3; d++)
            {
               TargetRange_Min[d] = MPI_Rank_X[d]*NX0[d]*amr->scale[0];
               TargetRange_Max[d] = TargetRange_Min[d] + NX0[d]*amr->scale[0];
            }
#           endif


            File = fopen( FileName, "rb" );
            fseek( File, Offset, SEEK_SET );

            for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
            {
//             d2. load the corner and son of this patch
               fread( Load_Cr_and_Son, sizeof(int), 4, File );

               for (int d=0; d<3; d++)    LoadCorner[d] *= rescale;


//             verify that the loaded patch is within the target range
#              ifdef LOAD_BALANCE
               if (  MPI_Rank == LB_Index2Rank( lv, LB_Corner2Index(lv,LoadCorner,CHECK_ON), CHECK_ON )  )
#              else
               if (  LoadCorner[0] >= TargetRange_Min[0]  &&  LoadCorner[0] < TargetRange_Max[0]  &&
                     LoadCorner[1] >= TargetRange_Min[1]  &&  LoadCorner[1] < TargetRange_Max[1]  &&
                     LoadCorner[2] >= TargetRange_Min[2]  &&  LoadCorner[2] < TargetRange_Max[2]     )
#              endif
               {
                  amr->pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2], -1, true, true, true );

//                d3. load the physical data if it is a leaf patch
                  if ( *LoadSon == -1 )
                  {
                     PID = amr->num[lv] - 1;

//                   d3-0. load the particle information (for leaf patches only)
#                    ifdef PARTICLE
                     fread( Load_NPar_and_GParID, sizeof(long), 2, File );

//                   note that we temporarily store GParID in the LB_Idx of Sg=1
//                   (since it's the only variable declared as long and it's useless anyway for Sg=1)
                     amr->patch[0][lv][PID]->NPar   = *Load_NPar;
                     amr->patch[1][lv][PID]->LB_Idx = *Load_GParID;

                     NParThisRank     += *Load_NPar;
                     MaxNParInOnePatch = MAX( MaxNParInOnePatch, *Load_NPar );
#                    endif

//                   d3-1. load the fluid variables
                     fread( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid, sizeof(real), CUBE(PS1)*NCOMP_TOTAL, File );

//                   d3-2. skip gravitational potential
#                    ifdef GRAVITY
                     if ( LoadPot )       fseek( File, CUBE(PS1)*sizeof(real), SEEK_CUR );
#                    endif

//                   d3-3. skip particle density
#                    ifdef PARTICLE
                     if ( LoadParDens )   fseek( File, CUBE(PS1)*sizeof(real), SEEK_CUR );
#                    endif

//                   d3-4. load magnetic field
#                    ifdef MHD
//                   skip the cell-centered data
                     if ( LoadCCMag )     fseek( File, CUBE(PS1)*NCOMP_MAG*sizeof(real), SEEK_CUR );

//                   load the face-centered data
                     fread( amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic, sizeof(real), PS1P1*SQR(PS1)*NCOMP_MAG, File );
#                    endif
                  } // if ( *LoadSon == -1 )
               } // within the target range

//             for the case that the patch is NOT within the target range
               else if ( *LoadSon == -1 )
               {
//                for particles, skip NPar and GParID as well
#                 ifdef PARTICLE
                  fseek( File, PatchDataSize+2*sizeof(long), SEEK_CUR );
#                 else
                  fseek( File, PatchDataSize, SEEK_CUR );
#                 endif
               }
            } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

            fclose( File );


//          d4. record the number of the real patches and the LB_IdxList_real
            for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = amr->num[lv];

#           ifdef LOAD_BALANCE
            if ( amr->LB->IdxList_Real         [lv] != NULL )   delete [] amr->LB->IdxList_Real         [lv];
            if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )   delete [] amr->LB->IdxList_Real_IdxTable[lv];

            amr->LB->IdxList_Real         [lv] = new long [ amr->NPatchComma[lv][1] ];
            amr->LB->IdxList_Real_IdxTable[lv] = new int  [ amr->NPatchComma[lv][1] ];

            for (int RPID=0; RPID<amr->NPatchComma[lv][1]; RPID++)
               amr->LB->IdxList_Real[lv][RPID] = amr->patch[0][lv][RPID]->LB_Idx;

            Mis_Heapsort( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );
#           endif // #ifdef LOAD_BALANCE

            Offset += DataSize[lv];

         } // if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )

         MPI_Barrier( MPI_COMM_WORLD );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

      } // for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)
   } // for (int lv=0; lv<NLv_Restart; lv++)


// get the total number of real patches at all ranks
   for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );



// e. load particles
// =================================================================================================
#  ifdef PARTICLE
   const long ParDataSize1v = amr->Par->NPar_Active_AllRank*sizeof(real);

   long  *NewParList = new long [MaxNParInOnePatch];
   real **ParBuf     = NULL;

   real NewParAtt[PAR_NATT_TOTAL];
   long GParID;
   int  NParThisPatch;

// be careful about using ParBuf returned from Aux_AllocateArray2D, which is set to NULL if MaxNParInOnePatch == 0
// --> for example, accessing ParBuf[0...PAR_NATT_STORED-1] will be illegal when MaxNParInOnePatch == 0
   Aux_AllocateArray2D( ParBuf, PAR_NATT_STORED, MaxNParInOnePatch );


// all particles are assumed to be synchronized with the base level
   NewParAtt[PAR_TIME] = Time[0];


// allocate particle repository
   amr->Par->InitRepo( NParThisRank, MPI_NRank );


// reset the total number of particles to be zero
// --> so particle repository is pre-allocated, but it contains no active particle yet
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;


// begin to load data
#  ifdef DEBUG_PARTICLE
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
#  endif

   for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "   Loading particle data, MPI ranks %4d -- %4d ... ",
                      TRanks, MIN(TRanks+RESTART_LOAD_NRANK-1, MPI_NRank-1) );

      if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )
      {
         File = fopen( FileName, "rb" );

         for (int lv=0; lv<NLEVEL; lv++)
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
//          note that we temporarily store GParID in the LBIdx of Sg=1
            NParThisPatch = amr->patch[0][lv][PID]->NPar;
            GParID        = amr->patch[1][lv][PID]->LB_Idx;

            if ( NParThisPatch > 0 )
            {
//             reset NPar=0 since we will call patch->AddParticle to update it again
               amr->patch[0][lv][PID]->NPar = 0;

//             load one particle attribute at a time
               for (int v=0; v<PAR_NATT_STORED; v++)
               {
                  fseek( File, FileOffset_Particle + v*ParDataSize1v + GParID*sizeof(real), SEEK_SET );

//                using ParBuf[v] here is safe since it's NOT called when NParThisPatch == 0
                  fread( ParBuf[v], sizeof(real), NParThisPatch, File );
               }

//             store particles to the particle repository (one particle at a time)
               for (int p=0; p<NParThisPatch; p++ )
               {
//                skip the last PAR_NATT_UNSTORED attributes since we do not store them on disk
                  for (int v=0; v<PAR_NATT_STORED; v++)  NewParAtt[v] = ParBuf[v][p];

                  NewParList[p] = amr->Par->AddOneParticle( NewParAtt );

#                 ifdef DEBUG_PARTICLE
                  if ( NewParList[p] >= NParThisRank )
                     Aux_Error( ERROR_INFO, "New particle ID (%ld) >= maximum allowed value (%ld) !!\n",
                                NewParList[p], NParThisRank );
#                 endif
               } // for (int p=0; p<NParThisPatch )

//             link particles to this patch
               const real *PType = amr->Par->Type;
#              ifdef DEBUG_PARTICLE
               char Comment[100];
               sprintf( Comment, "%s, PID %d, NPar %d", __FUNCTION__, PID, NParThisPatch );
               amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                                    PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#              else
               amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParList, &amr->Par->NPar_Lv[lv],
                                                    PType );
#              endif
            } // if ( amr->patch[0][lv][PID]->NPar > 0 )
         } // for PID, lv

         fclose( File );

         if ( amr->Par->NPar_AcPlusInac != NParThisRank )
            Aux_Error( ERROR_INFO, "total number of particles in the repository (%ld) != expect (%ld) !!\n",
                       amr->Par->NPar_AcPlusInac, NParThisRank );
      } // if ( MPI_Rank >= TRanks  &&  MPI_Rank < TRanks+RESTART_LOAD_NRANK )

      MPI_Barrier( MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int TRanks=0; TRanks<MPI_NRank; TRanks+=RESTART_LOAD_NRANK)


// free memory
   delete [] NewParList;
   Aux_DeallocateArray2D( ParBuf );
#  endif // #ifdef PARTICLE



// f-1. improve load balance
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


// fill up the data of non-leaf patches
// --> only necessary when restarting from a C-binary snapshot since it does not store non-leaf data
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT, _TOTAL, _MAG );

      LB_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT, _TOTAL, _MAG, NULL_INT );

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );
   }



// f-2. complete all levels for the case without LOAD_BALANCE
// ===================================================================================================================
#  else // #ifdef LOAD_BALANCE

   for (int lv=0; lv<NLEVEL; lv++)
   {
//    construct the relation "father <-> son" for the in-core computing
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
   } // for (int lv=0; lv<NLEVEL; lv++)


// fill up the data for top-level buffer patches
   Buf_GetBufferData( NLEVEL-1, amr->FluSg[NLEVEL-1], amr->MagSg[NLEVEL-1], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );


// fill up the data for patches that are not leaf patches
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
//    data restriction: lv+1 --> lv
      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT, _TOTAL, _MAG );

//    fill up the data in the buffer patches
      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );
   } // for (int lv=NLEVEL-2; lv>=0; lv--)

#  endif // #ifdef LOAD_BALANCE ... else ...


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByRestart



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_After_2000
// Description :  Load all simulation parameters from the RESTART file with format version >= 2000
//
// Note        :  All floating-point variables are declared as double after version 2000
//
// Parameter   :  File           : RESTART file pointer
//                FormatVersion  : Format version of the RESTART file
//                NLv_Restart    : NLEVEL recorded in the RESTART file
//                LoadPot        : Whether or not the RESTART file stores the potential data
//                LoadParDens    : Whether or not the RESTART file stores the particle (or total) density data
//                LoadCCMag      : Whether or not the RESTART file stores the cell-centered magnetic field
//                HeaderOffset_X : Offsets of different headers
//
// Return      :  NLv_Restart, LoadPot, LoadParDens, LoadCCMag (END_T and END_STEP may also be set to the original values)
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, int &NLv_Restart,
                                bool &LoadPot, bool &LoadParDens, bool &LoadCCMag,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, gamer_optimization, gamer_debug, timing, timing_solver;
   bool intel, float8, serial, overlap_mpi, openmp, store_pot_ghost, unsplit_gravity, particle;
   bool conserve_mass, laplacian_4th, self_interaction, laohu, support_hdf5, mhd;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch, ncomp_passive, gpu_arch;
   int  eos;

   fseek( File, HeaderOffset_Makefile, SEEK_SET );

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &gamer_optimization,         sizeof(bool),                    1,             File );
   fread( &gamer_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &gpu_arch,                   sizeof(int),                     1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );
   fread( &store_pot_ghost,            sizeof(bool),                    1,             File );
   fread( &unsplit_gravity,            sizeof(bool),                    1,             File );
   fread( &particle,                   sizeof(bool),                    1,             File );
   fread( &ncomp_passive,              sizeof(int),                     1,             File );
   fread( &conserve_mass,              sizeof(bool),                    1,             File );
   fread( &laplacian_4th,              sizeof(bool),                    1,             File );
   fread( &self_interaction,           sizeof(bool),                    1,             File );
   fread( &laohu,                      sizeof(bool),                    1,             File );
   fread( &support_hdf5,               sizeof(bool),                    1,             File );
   fread( &mhd,                        sizeof(bool),                    1,             File );
   fread( &eos,                        sizeof(int),                     1,             File );


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool   enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate_useless;
   bool   use_psolver_10to14;
   int    ncomp_fluid, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate;
   int    flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size;
   int    par_natt_stored, par_natt_user;
   double min_pres, max_error;

   fseek( File, HeaderOffset_Constant, SEEK_SET );

   fread( &ncomp_fluid,                sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_pres,                   sizeof(double),                  1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate_useless,      sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(double),                  1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size,             sizeof(int),                     1,             File );
   fread( &par_natt_stored,            sizeof(int),                     1,             File );
   fread( &par_natt_user,              sizeof(int),                     1,             File );


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   dummy_bool, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_user, opt__output_base, opt__output_pot;
   bool   opt__output_baseps, opt__timing_balance, opt__int_phase, opt__1st_flux_corr, opt__unit;
   bool   opt__output_cc_mag;
   int    nx0_tot[3], mpi_nrank, mpi_nrank_x[3], omp_nthread, regrid_count, opt__output_par_dens;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter_useless, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step, opt__1st_flux_corr_scheme;
   long   end_step;
   double lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, elbdm_planck_const, newton_g, sor_omega;
   double mg_tolerated_error, output_part_x, output_part_y, output_part_z, molecular_weight;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt, hubble0;
   double unit_l, unit_m, unit_t, unit_v, unit_d, unit_e, unit_p, unit_b;

   fseek( File, HeaderOffset_Parameter, SEEK_SET );

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &mpi_nrank,                  sizeof(int),                     1,             File );
   fread(  mpi_nrank_x,                sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &dummy_bool,                 sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(double),                  1,             File );
   fread( &gamma,                      sizeof(double),                  1,             File );
   fread( &minmod_coeff,               sizeof(double),                  1,             File );
   fread( &ep_coeff,                   sizeof(double),                  1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter_useless,   sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(double),                  1,             File );
   fread( &elbdm_planck_const,         sizeof(double),                  1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(double),                  1,             File );
   fread( &sor_omega,                  sizeof(double),                  1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(double),                  1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_user,           sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(double),                  1,             File );
   fread( &output_part_y,              sizeof(double),                  1,             File );
   fread( &output_part_z,              sizeof(double),                  1,             File );
   fread( &opt__timing_balance,        sizeof(bool),                    1,             File );
   fread( &opt__output_baseps,         sizeof(bool),                    1,             File );
   fread( &opt__1st_flux_corr,         sizeof(bool),                    1,             File );
   fread( &opt__1st_flux_corr_scheme,  sizeof(int),                     1,             File );
   fread( &opt__output_par_dens,       sizeof(int),                     1,             File );
   fread( &hubble0,                    sizeof(double),                  1,             File );
   fread( &opt__unit,                  sizeof(bool),                    1,             File );
   fread( &unit_l,                     sizeof(double),                  1,             File );
   fread( &unit_m,                     sizeof(double),                  1,             File );
   fread( &unit_t,                     sizeof(double),                  1,             File );
   fread( &unit_v,                     sizeof(double),                  1,             File );
   fread( &unit_d,                     sizeof(double),                  1,             File );
   fread( &unit_e,                     sizeof(double),                  1,             File );
   fread( &unit_p,                     sizeof(double),                  1,             File );
   fread( &molecular_weight,           sizeof(double),                  1,             File );
   fread( &opt__output_cc_mag,         sizeof(bool),                    1,             File );
   fread( &unit_b,                     sizeof(double),                  1,             File );


// set some default parameters
   if ( END_T < 0.0  &&  !OPT__RESTART_RESET )
   {
      END_T = end_t;
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %14.7e\n", "END_T", END_T );
   }

   if ( END_STEP < 0  &&  !OPT__RESTART_RESET )
   {
      END_STEP = end_step;
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "      NOTE : parameter %s is reset to %ld\n", "END_STEP", END_STEP );
   }

   if ( FormatVersion < 2210 )   opt__output_cc_mag = false;

   if ( FormatVersion < 2220 )   eos                = EOS_GAMMA;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "   Checking loaded parameters ...\n" );


      const bool Fatal    = true;
      const bool NonFatal = false;

//    d-1. check the simulation options and parameters defined in the Makefile
//    ========================================================================

//    errors
//    ------------------
      CompareVar( "MODEL", model, MODEL, Fatal );

#     ifdef GRAVITY
      if ( !gravity )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "GRAVITY", "OFF", "ON" );
#     else
      if (  gravity )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "GRAVITY", "ON", "OFF" );
#     endif

#     ifdef COMOVING
      if ( !comoving )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "COMOVING", "OFF", "ON" );
#     else
      if (  comoving )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "COMOVING", "ON", "OFF" );
#     endif

#     ifdef FLOAT8
      if ( !float8 )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
#     else
      if (  float8 )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
#     endif

#     ifdef PARTICLE
      if ( !particle )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "PARTICLE", "OFF", "ON" );
#     else
      if (  particle )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "PARTICLE", "ON", "OFF" );
#     endif

#     ifdef MHD
      if ( !mhd  ||  FormatVersion < 2210 )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "MHD", "OFF", "ON" );
#     else
      if (  mhd  &&  FormatVersion >= 2210 )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%s) != runtime (%s) !!\n", "MHD", "ON", "OFF" );
#     endif

      if ( nlevel > NLEVEL )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%d) > runtime (%d) (please set NLEVEL larger) !!\n",
                    "NLEVEL", nlevel, NLEVEL );

      CompareVar( "NCOMP_PASSIVE", ncomp_passive, NCOMP_PASSIVE, Fatal );



//    warnings
//    ------------------
      if ( OPT__DT_LEVEL != DT_LEVEL_SHARED  &&  !individual_timestep )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "INDIVIDUAL_TIMESTEP", "OFF", "ON" );
      else if ( OPT__DT_LEVEL == DT_LEVEL_SHARED  &&  individual_timestep )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "INDIVIDUAL_TIMESTEP", "ON", "OFF" );

#     ifdef GPU
      if ( !gpu )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "GPU", "OFF", "ON" );
#     else
      if (  gpu )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "GPU", "ON", "OFF" );
#     endif

#     ifdef GAMER_DEBUG
      if ( !gamer_debug )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "GAMER_DEBUG", "OFF", "ON" );
#     else
      if (  gamer_debug )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "GAMER_DEBUG", "ON", "OFF" );
#     endif

#     ifdef TIMING
      if ( !timing )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "TIMING", "OFF", "ON" );
#     else
      if (  timing )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "TIMING", "ON", "OFF" );
#     endif

#     ifdef TIMING_SOLVER
      if ( !timing_solver )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "TIMING_SOLVER", "OFF", "ON");
#     else
      if (  timing_solver )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "TIMING_SOLVER", "ON", "OFF");
#     endif

#     ifdef SERIAL
      if ( !serial )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "SERIAL", "OFF", "ON" );
#     else
      if (  serial )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "SERIAL", "ON", "OFF" );
#     endif

#     ifdef OVERLAP_MPI
      if ( !overlap_mpi )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "OVERLAP_MPI", "OFF", "ON" );
#     else
      if (  overlap_mpi )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "OVERLAP_MPI", "ON", "OFF" );
#     endif

#     ifdef OPENMP
      if ( !openmp )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "OPENMP", "OFF", "ON" );
#     else
      if (  openmp )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "OPENMP", "ON", "OFF" );
#     endif

#     ifdef GPU
      CompareVar( "GPU_ARCH", gpu_arch, GPU_ARCH, NonFatal );
#     endif

#     ifdef LOAD_BALANCE
      CompareVar( "LOAD_BALANCE", load_balance, LOAD_BALANCE, NonFatal );
#     else
      if ( load_balance )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%d) != runtime (%s) !!\n",
                      "LOAD_BALANCE", load_balance, "OFF" );
#     endif

#     ifdef SUPPORT_HDF5
      if ( !support_hdf5 )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "SUPPORT_HDF5", "OFF", "ON" );
#     else
      if (  support_hdf5 )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "SUPPORT_HDF5", "ON", "OFF" );
#     endif

#     ifdef LAOHU
      if ( !laohu )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "LAOHU", "OFF", "ON" );
#     else
      if (  laohu )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n", "LAOHU", "ON", "OFF" );
#     endif

      if ( nlevel < NLEVEL )
      {
         Aux_Message( stderr, "WARNING : %s : RESTART file (%d) < runtime (%d) !!\n", "NLEVEL", nlevel, NLEVEL );
         Aux_Message( stderr, "          --> Grid scale will be rescaled\n" );
      }

      CompareVar( "MAX_PATCH", max_patch, MAX_PATCH, NonFatal );

      CompareVar( "EOS",       eos,       EOS,       NonFatal );



//    d-2. check the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
//    ========================================================================
      CompareVar( "NCOMP_FLUID",             ncomp_fluid,            NCOMP_FLUID,                  Fatal );
      CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );

      CompareVar( "FLU_GHOST_SIZE",          flu_ghost_size,         FLU_GHOST_SIZE,            NonFatal );
      CompareVar( "FLU_BLOCK_SIZE_X",        flu_block_size_x,       FLU_BLOCK_SIZE_X,          NonFatal );
      CompareVar( "FLU_BLOCK_SIZE_Y",        flu_block_size_y,       FLU_BLOCK_SIZE_Y,          NonFatal );
#     if ( MODEL == HYDRO )
      CompareVar( "MIN_PRES",                min_pres,               MIN_PRES,                  NonFatal );
#     endif


//    check in GRAVITY
//    ----------------
#     ifdef GRAVITY

      CompareVar( "POT_SCHEME",              pot_scheme,             POT_SCHEME,                NonFatal );
      CompareVar( "POT_GHOST_SIZE",          pot_ghost_size,         POT_GHOST_SIZE,            NonFatal );
      CompareVar( "GRA_GHOST_SIZE",          gra_ghost_size,         GRA_GHOST_SIZE,            NonFatal );

#     ifdef POT_BLOCK_SIZE_X
      CompareVar( "POT_BLOCK_SIZE_X",        pot_block_size_x,       POT_BLOCK_SIZE_X,          NonFatal );
#     endif

#     ifdef POT_BLOCK_SIZE_Z
      CompareVar( "POT_BLOCK_SIZE_Z",        pot_block_size_z,       POT_BLOCK_SIZE_Z,          NonFatal );
#     endif

#     ifdef GRA_BLOCK_SIZE
      CompareVar( "GRA_BLOCK_SIZE",          gra_block_size,         GRA_BLOCK_SIZE,            NonFatal );
#     endif

#     if ( POT_SCHEME == SOR )
#     ifdef USE_PSOLVER_10TO14
      if ( !use_psolver_10to14 )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "USE_PSOLVER_10TO14", "OFF", "ON" );
#     else
      if (  use_psolver_10to14 )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "USE_PSOLVER_10TO14", "ON", "OFF" );
#     endif
#     endif // if ( POT_SCHEME == SOR )

#     ifdef STORE_POT_GHOST
      if ( !store_pot_ghost )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "STORE_POT_GHOST", "OFF", "ON" );
#     else
      if (  store_pot_ghost )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "STORE_POT_GHOST", "ON", "OFF" );
#     endif

#     ifdef UNSPLIT_GRAVITY
      if ( !unsplit_gravity )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "UNSPLIT_GRAVITY", "OFF", "ON" );
#     else
      if (  unsplit_gravity )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "UNSPLIT_GRAVITY", "ON", "OFF" );
#     endif

#     endif // #ifdef GRAVITY


//    check in HYDRO
//    ----------------
#     if ( MODEL == HYDRO )

      CompareVar( "FLU_SCHEME",              flu_scheme,             FLU_SCHEME,                NonFatal );

#     ifdef LR_SCHEME
      CompareVar( "LR_SCHEME",               lr_scheme,              LR_SCHEME,                 NonFatal );
#     endif

#     ifdef RSOLVER
      CompareVar( "RSOLVER",                 rsolver,                RSOLVER,                   NonFatal );
#     endif

#     ifdef CHECK_INTERMEDIATE
      CompareVar( "CHECK_INTERMEDIATE",      check_intermediate,     CHECK_INTERMEDIATE,        NonFatal );
#     endif

      CompareVar( "MAX_ERROR",               max_error,      (double)MAX_ERROR,                 NonFatal );

      if ( !enforce_positive )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "MIN_DENS/MIN_PRES", "OFF", "ON" );

#     ifdef CHAR_RECONSTRUCTION
      if ( !char_reconstruction )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "CHAR_RECONSTRUCTION", "OFF", "ON" );
#     else
      if (  char_reconstruction )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "CHAR_RECONSTRUCTION", "ON", "OFF" );
#     endif

#     ifdef HLL_NO_REF_STATE
      if ( !hll_no_ref_state )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "HLL_NO_REF_STATE", "OFF", "ON" );
#     else
      if (  hll_no_ref_state )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "HLL_NO_REF_STATE", "ON", "OFF" );
#     endif

#     ifdef HLL_INCLUDE_ALL_WAVES
      if ( !hll_include_all_waves )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "HLL_INCLUDE_ALL_WAVES", "OFF", "ON" );
#     else
      if (  hll_include_all_waves )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "HLL_INCLUDE_ALL_WAVES", "ON", "OFF" );
#     endif


//    check in ELBDM
//    ----------------
#     elif ( MODEL == ELBDM )
#     ifdef CONSERVE_MASS
      if ( !conserve_mass )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "CONSERVE_MASS", "OFF", "ON" );
#     else
      if (  conserve_mass )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "CONSERVE_MASS", "ON", "OFF" );
#     endif

#     ifdef LAPLACIAN_4TH
      if ( !laplacian_4th )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "LAPLACIAN_4TH", "OFF", "ON" );
#     else
      if (  laplacian_4th )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "LAPLACIAN_4TH", "ON", "OFF" );
#     endif

#     ifdef QUARTIC_SELF_INTERACTION
      if ( !self_interaction )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "QUARTIC_SELF_INTERACTION", "OFF", "ON" );
#     else
      if (  self_interaction )
         Aux_Message( stderr, "WARNING : %s : RESTART file (%s) != runtime (%s) !!\n",
                      "QUARTIC_SELF_INTERACTION", "ON", "OFF" );
#     endif

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL


//    check in PARTICLE
//    ------------------
#     ifdef PARTICLE
      if ( FormatVersion >= 2200 )
      CompareVar( "PAR_NATT_STORED",         par_natt_stored,              PAR_NATT_STORED,              Fatal );
      else
      CompareVar( "PAR_NATT_STORED",         par_natt_stored,              PAR_NATT_STORED,           NonFatal );

      if ( FormatVersion >= 2200 )
      CompareVar( "PAR_NATT_USER",           par_natt_user,                PAR_NATT_USER,                Fatal );
      else
      CompareVar( "PAR_NATT_USER",           par_natt_user,                PAR_NATT_USER,             NonFatal );

      if ( par_natt_user > 0  &&  FormatVersion < 2200  &&  MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : loading user-defined particle attributes (PAR_NATT_USER = %d) "
                              "from version < 2200 will likely fail !!\n", par_natt_user );
#     endif



//    d-3. check the simulation parameters recorded in the file "Input__Parameter"
//    =================================================================================================

//    errors
//    ------------------
      CompareVar( "BOX_SIZE",                box_size,                     BOX_SIZE,                     Fatal );
      CompareVar( "NX0_TOT[0]",              nx0_tot[0],                   NX0_TOT[0],                   Fatal );
      CompareVar( "NX1_TOT[1]",              nx0_tot[1],                   NX0_TOT[1],                   Fatal );
      CompareVar( "NX2_TOT[2]",              nx0_tot[2],                   NX0_TOT[2],                   Fatal );


//    warnings
//    ------------------
      CompareVar( "MPI_NRank",               mpi_nrank,                    MPI_NRank,                 NonFatal );
      CompareVar( "MPI_NRank_X[0]",          mpi_nrank_x[0],               MPI_NRank_X[0],            NonFatal );
      CompareVar( "MPI_NRank_X[1]",          mpi_nrank_x[1],               MPI_NRank_X[1],            NonFatal );
      CompareVar( "MPI_NRank_X[2]",          mpi_nrank_x[2],               MPI_NRank_X[2],            NonFatal );
      CompareVar( "OMP_NTHREAD",             omp_nthread,                  OMP_NTHREAD,               NonFatal );
      CompareVar( "END_T",                   end_t,                        END_T,                     NonFatal );
      CompareVar( "END_STEP",                end_step,                     END_STEP,                  NonFatal );
      CompareVar( "DT__FLUID",               dt__fluid,                    DT__FLUID,                 NonFatal );
      CompareVar( "OPT__DT_USER",            opt__dt_user,                 OPT__DT_USER,              NonFatal );
      CompareVar( "REGRID_COUNT",            regrid_count,                 REGRID_COUNT,              NonFatal );
      CompareVar( "FLAG_BUFFER_SIZE",        flag_buffer_size,             FLAG_BUFFER_SIZE,          NonFatal );
      CompareVar( "MAX_LEVEL",               max_level,                    MAX_LEVEL,                 NonFatal );
      CompareVar( "OPT__FLAG_RHO",           OPT__FLAG_RHO,                OPT__FLAG_RHO,             NonFatal );
      CompareVar( "OPT__FLAG_RHO_GRADIENT",  opt__flag_rho_gradient,       OPT__FLAG_RHO_GRADIENT,    NonFatal );
      CompareVar( "OPT__FLAG_USER",          opt__flag_user,               OPT__FLAG_USER,            NonFatal );
      CompareVar( "FLU_GPU_NPGROUP",         flu_gpu_npgroup,              FLU_GPU_NPGROUP,           NonFatal );
      CompareVar( "GPU_NSTREAM",             gpu_nstream,                  GPU_NSTREAM,               NonFatal );
      CompareVar( "OPT__FIXUP_FLUX",         opt__fixup_flux,              OPT__FIXUP_FLUX,           NonFatal );
      CompareVar( "OPT__FIXUP_RESTRICT",     opt__fixup_restrict,          OPT__FIXUP_RESTRICT,       NonFatal );
      CompareVar( "OPT__OVERLAP_MPI",        opt__overlap_mpi,             OPT__OVERLAP_MPI,          NonFatal );
      CompareVar( "OPT__INT_TIME",           opt__int_time,                OPT__INT_TIME,             NonFatal );
      CompareVar( "OPT__FLU_INT_SCHEME",     opt__flu_int_scheme,     (int)OPT__FLU_INT_SCHEME,       NonFatal );
      CompareVar( "OPT__REF_FLU_INT_SCHEME", opt__ref_flu_int_scheme, (int)OPT__REF_FLU_INT_SCHEME,   NonFatal );
      CompareVar( "OPT__OUTPUT_TOTAL",       opt__output_total,       (int)OPT__OUTPUT_TOTAL,         NonFatal );
      CompareVar( "OPT__OUTPUT_PART",        opt__output_part,        (int)OPT__OUTPUT_PART,          NonFatal );
      CompareVar( "OPT__OUTPUT_USER",        opt__output_user,             OPT__OUTPUT_USER,          NonFatal );
      CompareVar( "OPT__OUTPUT_BASEPS",      opt__output_baseps,           OPT__OUTPUT_BASEPS,        NonFatal );
      if ( OPT__OUTPUT_PART )
      CompareVar( "OPT__OUTPUT_BASE",        opt__output_base,             OPT__OUTPUT_BASE,          NonFatal );
#     ifdef PARTICLE
      if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PAR_TEXT ) {
#     else
      if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS ) {
#     endif
      CompareVar( "OPT__OUTPUT_MODE",        opt__output_mode,        (int)OPT__OUTPUT_MODE,          NonFatal );
      CompareVar( "OUTPUT_STEP",             output_step,                  OUTPUT_STEP,               NonFatal );
      CompareVar( "OUTPUT_DT",               output_dt,                    OUTPUT_DT,                 NonFatal );}
      if ( OPT__OUTPUT_PART ) {
      CompareVar( "OUTPUT_PART_X",           output_part_x,                OUTPUT_PART_X,             NonFatal );
      CompareVar( "OUTPUT_PART_Y",           output_part_y,                OUTPUT_PART_Y,             NonFatal );
      CompareVar( "OUTPUT_PART_Z",           output_part_z,                OUTPUT_PART_Z,             NonFatal );}
      CompareVar( "OPT__TIMING_BALANCE",     opt__timing_balance,          OPT__TIMING_BALANCE,       NonFatal );

#     ifdef COMOVING
      CompareVar( "OMEGA_M0",                omega_m0,                     OMEGA_M0,                  NonFatal );
      CompareVar( "DT__MAX_DELTA_A",         dt__max_delta_a,              DT__MAX_DELTA_A,           NonFatal );
      if ( FormatVersion >= 2111 )
      CompareVar( "HUBBLE0",                 hubble0,                      HUBBLE0,                   NonFatal );
      else if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : restart file does not have the parameter \"%s\" !!\n", "HUBBLE0" );
#     endif

#     ifdef LOAD_BALANCE
      CompareVar( "LB->WLI_Max",             lb_wli_max,                   amr->LB->WLI_Max,          NonFatal );
#     endif

#     ifdef GRAVITY
      CompareVar( "DT__GRAVITY",             dt__gravity,                  DT__GRAVITY,               NonFatal );
      CompareVar( "NEWTON_G",                newton_g,                     NEWTON_G,                  NonFatal );
      CompareVar( "POT_GPU_NPGROUP",         pot_gpu_npgroup,              POT_GPU_NPGROUP,           NonFatal );
      CompareVar( "OPT__GRA_P5_GRADIENT",    opt__gra_p5_gradient,         OPT__GRA_P5_GRADIENT,      NonFatal );
      CompareVar( "OPT__POT_INT_SCHEME",     opt__pot_int_scheme,     (int)OPT__POT_INT_SCHEME,       NonFatal );
      CompareVar( "OPT__RHO_INT_SCHEME",     opt__rho_int_scheme,     (int)OPT__RHO_INT_SCHEME,       NonFatal );
      CompareVar( "OPT__GRA_INT_SCHEME",     opt__gra_int_scheme,     (int)OPT__GRA_INT_SCHEME,       NonFatal );
      CompareVar( "OPT__REF_POT_INT_SCHEME", opt__ref_pot_int_scheme, (int)OPT__REF_POT_INT_SCHEME,   NonFatal );
      CompareVar( "OPT__OUTPUT_POT",         opt__output_pot,              OPT__OUTPUT_POT,           NonFatal );
#     if   ( POT_SCHEME == SOR )
      CompareVar( "SOR_OMEGA",               sor_omega,                    SOR_OMEGA,                 NonFatal );
      CompareVar( "SOR_MAX_ITER",            sor_max_iter,                 SOR_MAX_ITER,              NonFatal );
      CompareVar( "SOR_MIN_ITER",            sor_min_iter,                 SOR_MIN_ITER,              NonFatal );
#     elif ( POT_SCHEME == MG )
      CompareVar( "MG_MAX_ITER",             mg_max_iter,                  MG_MAX_ITER,               NonFatal );
      CompareVar( "MG_NPRE_SMOOTH",          mg_npre_smooth,               MG_NPRE_SMOOTH,            NonFatal );
      CompareVar( "MG_NPOST_SMOOTH",         mg_npost_smooth,              MG_NPOST_SMOOTH,           NonFatal );
      CompareVar( "MG_TOLERATED_ERROR",      mg_tolerated_error,           MG_TOLERATED_ERROR,        NonFatal );
#     endif // POT_SCHEME
#     endif // #ifdef GRAVITY

#     if   ( MODEL == HYDRO )
      CompareVar( "OPT__FLAG_PRES_GRADIENT", opt__flag_pres_gradient,      OPT__FLAG_PRES_GRADIENT,   NonFatal );
      CompareVar( "GAMMA",                   gamma,                        GAMMA,                     NonFatal );
      if ( FormatVersion >= 2111 )
      CompareVar( "MOLECULAR_WEIGHT",        molecular_weight,             MOLECULAR_WEIGHT,          NonFatal );
      else
      Aux_Message( stderr, "WARNING : restart file does not have the parameter \"%s\" !!\n", "MOLECULAR_WEIGHT" );
      CompareVar( "MINMOD_COEFF",            minmod_coeff,                 MINMOD_COEFF,              NonFatal );
      CompareVar( "OPT__LR_LIMITER",         opt__lr_limiter,         (int)OPT__LR_LIMITER,           NonFatal );

//    convert OPT__1ST_FLUX_CORR to bool to be consistent with the old format where OPT__1ST_FLUX_CORR is bool instead of int
      CompareVar( "OPT__1ST_FLUX_CORR",        opt__1st_flux_corr,        (bool)OPT__1ST_FLUX_CORR,        NonFatal );
      CompareVar( "OPT__1ST_FLUX_CORR_SCHEME", opt__1st_flux_corr_scheme, (int )OPT__1ST_FLUX_CORR_SCHEME, NonFatal );

#     elif ( MODEL == ELBDM )
      CompareVar( "DT__PHASE",               dt__phase,                    DT__PHASE,                 NonFatal );
      CompareVar( "OPT__FLAG_ENGY_DENSITY",  opt__flag_engy_density,       OPT__FLAG_ENGY_DENSITY,    NonFatal );
      CompareVar( "OPT__INT_PHASE",          opt__int_phase,               OPT__INT_PHASE,            NonFatal );
      CompareVar( "ELBDM_MASS",              elbdm_mass,                   ELBDM_MASS,                NonFatal );
      CompareVar( "ELBDM_PLANCK_CONST",      elbdm_planck_const,           ELBDM_PLANCK_CONST,        NonFatal );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif

#     ifdef PARTICLE
      CompareVar( "OPT__OUTPUT_PAR_DENS",    opt__output_par_dens,    (int)OPT__OUTPUT_PAR_DENS,      NonFatal );
#     endif

      if ( FormatVersion >= 2111 ) {
      CompareVar( "OPT__UNIT",               opt__unit,                    OPT__UNIT,                 NonFatal );
      CompareVar( "UNIT_L",                  unit_l,                       UNIT_L,                    NonFatal );
      CompareVar( "UNIT_M",                  unit_m,                       UNIT_M,                    NonFatal );
      CompareVar( "UNIT_T",                  unit_t,                       UNIT_T,                    NonFatal );
      CompareVar( "UNIT_V",                  unit_v,                       UNIT_V,                    NonFatal );
      CompareVar( "UNIT_D",                  unit_d,                       UNIT_D,                    NonFatal );
      CompareVar( "UNIT_E",                  unit_e,                       UNIT_E,                    NonFatal );
      CompareVar( "UNIT_P",                  unit_p,                       UNIT_P,                    NonFatal ); }
      else if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : restart file does not have any information about the code units !!\n" );

#     ifdef MHD
      CompareVar( "OPT__OUTPUT_CC_MAG",      opt__output_cc_mag,           OPT__OUTPUT_CC_MAG,        NonFatal );
      CompareVar( "UNIT_B",                  unit_b,                       UNIT_B,                    NonFatal );
#     endif

      Aux_Message( stdout, "   Checking loaded parameters ... done\n" );

   } // if ( MPI_Rank == 0 )


// set the returned variables
#  ifdef GRAVITY
   LoadPot     = opt__output_pot;
#  else
   LoadPot     = false;
#  endif
#  ifdef PARTICLE
   LoadParDens = ( opt__output_par_dens != (int)PAR_OUTPUT_DENS_NONE );
#  else
   LoadParDens = false;
#  endif
#  ifdef MHD
   LoadCCMag   = opt__output_cc_mag;
#  else
   LoadCCMag   = false;
#  endif
   NLv_Restart = nlevel;

} // FUNCTION : Load_Parameter_After_2000



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareVar
// Description :  Compare the input variables
//
// Note        :  This function is overloaded to work with different data types
//
// Parameter   :  VarName     : Name of the target variable
//                RestartVar  : Variable loaded from the RESTART file
//                RuntimeVar  : Variable loaded from the Input__Parameter
//                Fatal       : Whether or not the difference between RestartVar and RuntimeVar is fatal
//                              --> true  : terminate the program if the input variables are different
//                                  false : display warning message if the input variables are different
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const bool RestartVar, const bool RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%d) != runtime (%d) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (bool)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "int"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%d) != runtime (%d) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (int)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "long"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const long RestartVar, const long RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%ld) != runtime (%ld) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%ld) != runtime (%ld) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (long)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "float"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const float RestartVar, const float RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (float)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "double"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
         Aux_Error( ERROR_INFO, "%s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                    VarName, RestartVar, RuntimeVar );
      else
         Aux_Message( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                      VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (double)
