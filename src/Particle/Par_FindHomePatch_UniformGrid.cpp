#include "GAMER.h"

#ifdef PARTICLE

static long ParPos2LBIdx( const int lv, const real ParPos[] );
static void SendParticle2HomeRank( const int lv, const bool OldParOnly,
                                   const long NNewPar, real *NewParAtt[PAR_NATT_TOTAL] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_FindHomePatch_UniformGrid
// Description :  Find the home patches of all particles on a fully occupied level
//
// Note        :  1. The target level must be fully occupied by patches
//                   --> So either "lv=0" or "lv-1 is fully refined"
//                2. Invoked by Init_UniformGrid(), Init_BaseLevel(), and Par_AddParticleAfterInit()
//                3. After calling this function, the amr->Par structure will be reconstructed and
//                   all particles will be associated with their home patches on lv
//
// Parameter   :  lv         : Target level
//                OldParOnly : true  --> only redistribute particles already exist in the current repository
//                                       --> Ignore NNewPar and NewParAtt
//                             false --> only redistribute newly added particles specified by NNewPar and NewParAtt
//                                       --> Particles already exist in the current repository will not be redistributed
//                NNewPar    : Number of new particles to be added                       (for OldParOnly==false only)
//                NewParAtt  : Pointer array storing the data of new particle attributes (for OldParOnly==false only)
//                             --> Format: real *NewParAtt[PAR_NATT_TOTAL]
//                             --> Must be deallocated manually after invoking this function
//
// Return      :  1. amr->Par
//                2. NPar, ParListSize, and ParList[] of all real patches on lv
//-------------------------------------------------------------------------------------------------------
void Par_FindHomePatch_UniformGrid( const int lv, const bool OldParOnly,
                                    const long NNewPar, real *NewParAtt[PAR_NATT_TOTAL] )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );

   if ( !OldParOnly )
   {
      if ( NNewPar < 0 )
         Aux_Error( ERROR_INFO, "NNewPar (%ld) < 0 !!\n", NNewPar );

      if ( NNewPar > 0  &&  NewParAtt[PAR_MASS] == NULL )
         Aux_Error( ERROR_INFO, "NewParAtt[PAR_MASS] == NULL !!\n" );
   }
#  endif


// record the number of old particles before it is overwritten by SendParticle2HomeRank()
   const long NOldPar = ( OldParOnly ) ? 0 : amr->Par->NPar_AcPlusInac;


// 1. redistribute all particles to their home ranks
   SendParticle2HomeRank( lv, OldParOnly, NNewPar, NewParAtt );


// 2. find the home patch
   const long  NewParID0 = NOldPar;
   const long  NTarPar   = amr->Par->NPar_AcPlusInac - NOldPar;
   const real *Pos[3]    = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *PType     = amr->Par->Type;
   const int   NReal     = amr->NPatchComma[lv][1];

   real TParPos[3];

   long *HomeLBIdx          = new long [NTarPar];
   int  *HomeLBIdx_IdxTable = new int  [NTarPar];
   int  *HomePID            = new int  [NTarPar];
   int  *MatchIdx           = new int  [NTarPar];

// get the load-balance index of the particle's home patch
   for (long t=0; t<NTarPar; t++)
   {
      for (int d=0; d<3; d++)    TParPos[d] = Pos[d][ NewParID0 + t ];
      HomeLBIdx[t] = ParPos2LBIdx( lv, TParPos );
   }

// sort the LBIdx list
   Mis_Heapsort( NTarPar, HomeLBIdx, HomeLBIdx_IdxTable );

// construct and sort the LBIdx list of all real patches
// --> do not use amr->LB->IdxList_Real[] since it may not be constructed yet
   long *RealPatchLBIdx          = new long [NReal];
   int  *RealPatchLBIdx_IdxTable = new int  [NReal];

   for (int PID=0; PID<NReal; PID++)   RealPatchLBIdx[PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NReal, RealPatchLBIdx, RealPatchLBIdx_IdxTable );

// LBIdx --> home (real) patch indices
   Mis_Matching_int( NReal, RealPatchLBIdx, NTarPar, HomeLBIdx, MatchIdx );

// check: every particle must have a home patch
#  ifdef DEBUG_PARTICLE
   for (long t=0; t<NTarPar; t++)
   {
      if ( MatchIdx[t] == -1 )
      {
         const long ParID = NewParID0 + HomeLBIdx_IdxTable[t];

         Aux_Error( ERROR_INFO, "lv %d, ParID %ld, ParPos (%14.7e, %14.7e, %14.7e) --> found no home patch !!\n",
                    lv, ParID, Pos[0][ParID], Pos[1][ParID], Pos[2][ParID] );
      }
   }
#  endif

// record the home patch indices
   for (long t=0; t<NTarPar; t++)
   {
      const long t0  = HomeLBIdx_IdxTable[t];
      const int  PID = RealPatchLBIdx_IdxTable[ MatchIdx[t] ];

      HomePID[t0] = PID;
   }


// 3. count the number of particles in each patch and allocate the particle list
   if ( OldParOnly )
   for (int PID=0; PID<NReal; PID++)   amr->patch[0][lv][PID]->ParListSize = 0;

   for (long t=0; t<NTarPar; t++)
      amr->patch[0][lv][ HomePID[t] ]->ParListSize ++;

   for (int PID=0; PID<NReal; PID++)
   {
      if ( OldParOnly )
      {
         free( amr->patch[0][lv][PID]->ParList );
         amr->patch[0][lv][PID]->ParList = NULL;
      }

      if ( amr->patch[0][lv][PID]->ParListSize > 0 )
      {
//       use realloc to preserve the old particle list for OldParOnly
         amr->patch[0][lv][PID]->ParList = (long*)realloc( amr->patch[0][lv][PID]->ParList,
                                                           amr->patch[0][lv][PID]->ParListSize*sizeof(long) );
      }
   }


// 4. associate particles with their home patches
   if ( OldParOnly )    amr->Par->NPar_Lv[lv] = 0;

// no OpenMP since AddParticle() will modify amr->Par->NPar_Lv[]
   for (long t=0; t<NTarPar; t++)
   {
      const long ParID = NewParID0 + t;
      const int  PID   = HomePID[t];

#     ifdef DEBUG_PARTICLE
      amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv], Ptype,
                                           Pos, amr->Par->NPar_AcPlusInac, __FUNCTION__ );
#     else
      amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv], PType);
#     endif
   }


   delete [] HomeLBIdx;
   delete [] HomeLBIdx_IdxTable;
   delete [] HomePID;
   delete [] MatchIdx;
   delete [] RealPatchLBIdx;
   delete [] RealPatchLBIdx_IdxTable;

} // FUNCTION : Par_FindHomePatch_UniformGrid



//-------------------------------------------------------------------------------------------------------
// Function    :  SendParticle2HomeRank
// Description :  Redistribute all particles to their home ranks
//
// Note        :  1. amr->LB->CutPoint[lv][] must be set in advance (e.g., by calling LB_SetCutPoint())
//                2. Invoked by Par_FindHomePatch_UniformGrid()
//                3. Inactive particles will NOT be redistributed
//                4. Currently only redistribute the Time, Type, Mass, Pos, Vel, and all passive variables
//                   of particles
//                   --> Acc is not redistributed since it has not been initialized when calling
//                       Par_FindHomePatch_UniformGrid()
//
// Parameter   :  lv         : Target level
//                OldParOnly : true  --> only redistribute particles already exist in the current repository
//                                       --> Ignore NNewPar and NewParAtt
//                             false --> only redistribute newly added particles specified by NNewPar and NewParAtt
//                                       --> Particles already exist in the current repository will not be redistributed
//                NNewPar    : Number of new particles to be added                       (for OldParOnly==false only)
//                NewParAtt  : Pointer array storing the data of new particle attributes (for OldParOnly==false only)
//                             --> Format: real *NewParAtt[PAR_NATT_TOTAL]
//                             --> Must be deallocated manually after invoking this function
//
// Return      :  1. amr->Par
//                2. NPar, ParListSize and ParList[] of all real patches on lv
//-------------------------------------------------------------------------------------------------------
void SendParticle2HomeRank( const int lv, const bool OldParOnly,
                            const long NNewPar, real *NewParAtt[PAR_NATT_TOTAL] )
{

   real *Mass   = NULL;
   real *Pos[3] = { NULL, NULL, NULL };

   int Send_Count[MPI_NRank], Recv_Count[MPI_NRank], Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];
   int Send_Count_Sum, Recv_Count_Sum;

   if ( OldParOnly )
   {
      Mass   = amr->Par->Mass;
      Pos[0] = amr->Par->PosX;
      Pos[1] = amr->Par->PosY;
      Pos[2] = amr->Par->PosZ;
   }

   else
   {
      Mass   = NewParAtt[PAR_MASS];
      Pos[0] = NewParAtt[PAR_POSX];
      Pos[1] = NewParAtt[PAR_POSY];
      Pos[2] = NewParAtt[PAR_POSZ];
   }

   for (int r=0; r<MPI_NRank; r++)  Send_Count[r] = 0;


// 1. get the target MPI rank of each particle
   const long NTarPar = ( OldParOnly ) ? amr->Par->NPar_AcPlusInac : NNewPar;

   int *TRank = new int [NTarPar];
   real TParPos[3];
   long LBIdx;

   for (long ParID=0; ParID<NTarPar; ParID++)
   {
//    TRank is set to -1 for inactive particles
      if ( Mass[ParID] < (real)0.0 )
      {
         TRank[ParID] = -1;
         continue;
      }

//    get the load-balance index of the particle's home patch
      for (int d=0; d<3; d++)    TParPos[d] = Pos[d][ParID];
      LBIdx = ParPos2LBIdx( lv, TParPos );

//    record the home patch
#     ifdef SERIAL
      TRank[ParID] = 0;
#     else
      TRank[ParID] = LB_Index2Rank( lv, LBIdx, CHECK_ON );
#     endif
      Send_Count[ TRank[ParID] ] ++;
   } // for (long ParID=0; ParID<NTarPar; ParID++)


// 2. construct the MPI send and recv data list
   MPI_Alltoall( Send_Count, 1, MPI_INT, Recv_Count, 1, MPI_INT, MPI_COMM_WORLD );

   Send_Disp[0] = 0;
   Recv_Disp[0] = 0;

   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + Send_Count[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + Recv_Count[r-1];
   }

   Send_Count_Sum = Send_Disp[ MPI_NRank-1 ] + Send_Count[ MPI_NRank-1 ];
   Recv_Count_Sum = Recv_Disp[ MPI_NRank-1 ] + Recv_Count[ MPI_NRank-1 ];

#  ifdef DEBUG_PARTICLE
   if ( OldParOnly  &&  Send_Count_Sum != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "Send_Count_Sum (%d) != NPar_Active (%ld) !!\n", Send_Count_Sum, amr->Par->NPar_Active );
#  endif


// 3. redistribute particle attributes (one attribute at a time to save memory)
   const long NOldPar            = ( OldParOnly ) ?       NULL_INT : amr->Par->NPar_AcPlusInac;
   const long UpdatedParListSize = ( OldParOnly ) ? Recv_Count_Sum : amr->Par->NPar_AcPlusInac + Recv_Count_Sum;

   int Offset[MPI_NRank];

   real *SendBuf = new real [Send_Count_Sum];
   real *RecvBuf = NULL;

// 3-1. record attribute pointers
   real *SendAttPtr[PAR_NATT_TOTAL], **OldAttPtrPtr[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
      SendAttPtr  [v] = ( OldParOnly ) ? amr->Par->Attribute[v] : NewParAtt[v];
      OldAttPtrPtr[v] = &amr->Par->Attribute[v];
   }

   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
//    3-2. initialize the array offsets of all target ranks
      for (int r=0; r<MPI_NRank; r++)  Offset[r] = Send_Disp[r];

//    3-3. prepare send buffer (skip inactive particles)
      for (long ParID=0; ParID<NTarPar; ParID++)
         if ( TRank[ParID] != -1 )  SendBuf[ Offset[TRank[ParID]] ++ ] = SendAttPtr[v][ParID];

//    3-4. free/allocate the old/new particle arrays and set the recv buffer
      if ( OldParOnly )
      {
         free( SendAttPtr[v] );

         *(OldAttPtrPtr[v]) = (real*)malloc( UpdatedParListSize*sizeof(real) );
         RecvBuf            = *(OldAttPtrPtr[v]);
      }

      else
      {
         *(OldAttPtrPtr[v]) = (real*)realloc( *(OldAttPtrPtr[v]), UpdatedParListSize*sizeof(real) );
         RecvBuf            = *(OldAttPtrPtr[v]) + NOldPar;
      }

//    3-5. redistribute data
#     ifdef FLOAT8
      MPI_Alltoallv( SendBuf, Send_Count, Send_Disp, MPI_DOUBLE, RecvBuf, Recv_Count, Recv_Disp, MPI_DOUBLE, MPI_COMM_WORLD );
#     else
      MPI_Alltoallv( SendBuf, Send_Count, Send_Disp, MPI_FLOAT,  RecvBuf, Recv_Count, Recv_Disp, MPI_FLOAT,  MPI_COMM_WORLD );
#     endif
   } // for (int v=0; v<PAR_NATT_TOTAL; v++)


// 4. reset particle parameters
   if ( OldParOnly )
   {
      free( amr->Par->InactiveParList );

      amr->Par->NPar_AcPlusInac     = (long)Recv_Count_Sum;
      amr->Par->NPar_Active         = (long)Recv_Count_Sum;
      amr->Par->NPar_Inactive       = 0;                       // since we don't redistribute inactive particles
      amr->Par->ParListSize         = UpdatedParListSize;
      amr->Par->InactiveParListSize = (long)MAX( 1, UpdatedParListSize/100 );
      amr->Par->InactiveParList     = (long*)malloc( amr->Par->InactiveParListSize*sizeof(long) );
   }

   else
   {
      amr->Par->NPar_AcPlusInac    += (long)Recv_Count_Sum;
      amr->Par->NPar_Active        += (long)Recv_Count_Sum;
      amr->Par->ParListSize         = UpdatedParListSize;

//    update the total number of active particles in all MPI ranks
      MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   } // if ( OldParOnly ) ... else ...


// 5. reset attribute pointers
   amr->Par->Mass = amr->Par->Attribute[PAR_MASS];
   amr->Par->PosX = amr->Par->Attribute[PAR_POSX];
   amr->Par->PosY = amr->Par->Attribute[PAR_POSY];
   amr->Par->PosZ = amr->Par->Attribute[PAR_POSZ];
   amr->Par->VelX = amr->Par->Attribute[PAR_VELX];
   amr->Par->VelY = amr->Par->Attribute[PAR_VELY];
   amr->Par->VelZ = amr->Par->Attribute[PAR_VELZ];
   amr->Par->Time = amr->Par->Attribute[PAR_TIME];
   amr->Par->Type = amr->Par->Attribute[PAR_TYPE];
#  ifdef STORE_PAR_ACC
   amr->Par->AccX = amr->Par->Attribute[PAR_ACCX];
   amr->Par->AccY = amr->Par->Attribute[PAR_ACCY];
   amr->Par->AccZ = amr->Par->Attribute[PAR_ACCZ];
#  endif


// 6. check
#  ifdef DEBUG_PARTICLE
   const long NewParID0 = ( OldParOnly ) ? 0 : NOldPar;

   for (long ParID=NewParID0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    there should be no inactive particles
      if ( amr->Par->Mass[ParID] < 0.0 )
         Aux_Error( ERROR_INFO, "ParMass[%ld] = %14.7e < 0.0 !!\n", ParID, amr->Par->Mass[ParID] );
   }
#  endif


// free memory
   delete [] TRank;
   delete [] SendBuf;

} // FUNCTION : SendParticle2HomeRank



//-------------------------------------------------------------------------------------------------------
// Function    :  ParPos2LBIdx
// Description :  Convert the input particle position to the load-balance index (LBIdx) of its home patch
//
// Note        :  1. Input particle position must lie within the simulation box
//                   --> Does not check periodicity
//                2. Invoked by SendParticle2HomeRank() and Par_FindHomePatch_UniformGrid()
//
// Parameter   :  lv     : Target level
//                ParPos : Particle position
//
// Return      :  Load-balance index (LBIdx)
//-------------------------------------------------------------------------------------------------------
long ParPos2LBIdx( const int lv, const real ParPos[] )
{

// check
#  ifdef DEBUG_PARTICLE
   for (int d=0; d<3; d++)
   {
      if ( ParPos[d] < amr->BoxEdgeL[d]  ||  ParPos[d] >= amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "ParPos[%d] = %14.7e lies outside the simulation box (BoxEdgeL/R = %14.7e/%14.7e) !!\n",
                    d, ParPos[d], amr->BoxEdgeL[d], amr->BoxEdgeR[d] );
   }
#  endif


// calculate the home patch corner
   const double dh_min        = amr->dh[TOP_LEVEL];
   const double _PatchPhySize = 1.0/( PS1*amr->dh[lv] );
   const int    PatchScale    = PS1*amr->scale[lv];

   double PatchEdgeL, PatchEdgeR;
   int    Cr[3];

   for (int d=0; d<3; d++)
   {
      Cr[d] = (int)floor(  ( (double)ParPos[d] - amr->BoxEdgeL[d] )*_PatchPhySize  )*PatchScale;

//    check the home patch corner carefully to prevent from any issue resulting from round-off errors
//    --> make sure to adopt the same procedure of calculating the patch left/right edges as Patch.h
      PatchEdgeL = (double)( Cr[d]              )*dh_min;
      PatchEdgeR = (double)( Cr[d] + PatchScale )*dh_min;

      if      ( ParPos[d] <  PatchEdgeL )    Cr[d] -= PatchScale;
      else if ( ParPos[d] >= PatchEdgeR )    Cr[d] += PatchScale;

//    check if the target home patch does enclose the given particle position
#     ifdef DEBUG_PARTICLE
      PatchEdgeL = (double)( Cr[d]              )*dh_min;
      PatchEdgeR = (double)( Cr[d] + PatchScale )*dh_min;

      if ( ParPos[d] < PatchEdgeL  ||  ParPos[d] >= PatchEdgeR )
         Aux_Error( ERROR_INFO, "incorrect home patch (dim %d, ParPos %14.7e, PatchEdgeL/R %14.7e/%14.7e) !!\n",
                    d, ParPos[d], PatchEdgeL, PatchEdgeR );

      if ( PatchEdgeL < amr->BoxEdgeL[d]  ||  PatchEdgeR > amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "incorrect home patch (dim %d, PatchEdgeL/R %14.7e/%14.7e, BoxEdgeL/R %14.7e/%14.7e) !!\n",
                    d, PatchEdgeL, PatchEdgeR, amr->BoxEdgeL[d], amr->BoxEdgeR[d] );
#     endif
   } // for (int d=0; d<3; d++)


// corner --> LBIdx and return
   const long LBIdx = LB_Corner2Index( lv, Cr, CHECK_ON );

   return LBIdx;

} // FUNCTION : ParPos2LBIdx



#endif // #ifdef PARTICLE
