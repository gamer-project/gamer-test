#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_ExchangeParticleBetweenPatch
// Description :  Exchange particles between patches in different MPI ranks
//
// Note        :  1. Information of both the send patches "Send_NPatchTotal, Send_PIDList, Send_NPatchEachRank"
//                   and the corresponding recv patches "Recv_NPatchTotal, Recv_PIDList, Recv_NPatchEachRank" must be
//                   provided. The mapping between send and recv patches can be calculated in advance by
//                   Par_LB_MapBuffer2RealPatch() if send patches are buffer and recv patches are real
//                2. All Target patches (those in Send_PIDList[] and Recv_PIDList[]) must be patches at the same level "lv"
//                3. Currently this function exchanges ALL particle attributes
//                   --> But it can be generalized to work with arbitrary particle attributes
//                4. This function is called by Par_PassParticle2Sibling(), Par_PassParticle2Son_MultiPatch(),
//                   and Par_LB_Refine_SendParticle2Father()
//
// Parameter   :  lv                  : Target refinement level
//                Send_NPatchTotal    : Total number of patches in Send_PIDList
//                Send_PIDList        : Patch indices to send particles
//                Send_NPatchEachRank : Number of patches to send particles to each rank
//                Recv_NPatchTotal    : Total number of patches in Recv_PIDList
//                Recv_PIDList        : Patch indices to receive particles
//                Recv_NPatchEachRank : Number of patches to receive particles from each rank
//                Timer               : Timer used by Par_LB_SendParticleData()
//                Timer_Comment       : String used by Par_LB_SendParticleData()
//
// Return      :  New particles will be added to the particle repository of this rank and linked to the
//                target recv patches
//-------------------------------------------------------------------------------------------------------
void Par_LB_ExchangeParticleBetweenPatch( const int lv,
                                          const int Send_NPatchTotal, const int *Send_PIDList, int *Send_NPatchEachRank,
                                          const int Recv_NPatchTotal, const int *Recv_PIDList, int *Recv_NPatchEachRank,
                                          Timer_t *Timer, const char *Timer_Comment )
{

// nothing to do for levels above MAX_LEVEL
   if ( lv > MAX_LEVEL )  return;


// check
#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv >= NLEVEL )     Aux_Error( ERROR_INFO, "incorrect target level (%d) !!\n", lv );

   if      ( Send_NPatchTotal < 0 )    Aux_Error( ERROR_INFO, "Send_NPatchTotal = %d < 0 !!\n", Send_NPatchTotal );
   else if ( Send_NPatchTotal > 0 )
   {
      if ( Send_PIDList == NULL )
         Aux_Error( ERROR_INFO, "Send_PIDList == NULL (Send_NPatchTotal = %d) !!\n", Send_NPatchTotal );

      if ( Send_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Send_NPatchEachRank == NULL (Send_NPatchTotal = %d) !!\n", Send_NPatchTotal );
   }

   if      ( Recv_NPatchTotal < 0 )    Aux_Error( ERROR_INFO, "Recv_NPatchTotal = %d < 0 !!\n", Recv_NPatchTotal );
   else if ( Recv_NPatchTotal > 0 )
   {
      if ( Recv_PIDList == NULL )
         Aux_Error( ERROR_INFO, "Recv_PIDList == NULL (Recv_NPatchTotal = %d) !!\n", Recv_NPatchTotal );

      if ( Recv_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Recv_NPatchEachRank == NULL (Recv_NPatchTotal = %d) !!\n", Recv_NPatchTotal );
   }

   const int *Check_PIDList    [2] = { Send_PIDList, Recv_PIDList };
   const int  Check_NPatchTotal[2] = { Send_NPatchTotal, Recv_NPatchTotal };

   for (int m=0; m<2; m++)
   for (int t=0; t<Check_NPatchTotal[m]; t++)
   {
      const int PID = Check_PIDList[m][t];

      if ( PID < 0  ||  PID >= amr->num[lv] )
         Aux_Error( ERROR_INFO, "This is NOT a correct patch index (lv %d, m %d, t %d, PID %d, NReal %d, NTotal %d) !!\n",
                    lv, m, t, PID, amr->NPatchComma[lv][1], amr->num[lv] );

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );

      if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar_Copy );

      for (int v=0; v<4; v++)
      if ( amr->patch[0][lv][PID]->ParMassPos_Copy[v] != NULL )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d, ParMassPos_Copy[%d] != NULL !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar_Copy, v );
   } // for m, t

// check Recv_NPatchEachRank and Recv_NPatchTotal
   int Recv_NPatchEachRank_Check[MPI_NRank], Recv_NPatchTotal_Check=0;

   MPI_Alltoall( Send_NPatchEachRank, 1, MPI_INT, Recv_NPatchEachRank_Check, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   if ( Recv_NPatchEachRank[r] != Recv_NPatchEachRank_Check[r] )
      Aux_Error( ERROR_INFO, "lv %d, Recv_NPatchEachRank[%d] (%d) != expected (%d) !!\n",
                 lv, r, Recv_NPatchEachRank[r], Recv_NPatchEachRank_Check[r] );

   for (int r=0; r<MPI_NRank; r++)  Recv_NPatchTotal_Check += Recv_NPatchEachRank_Check[r];

   if ( Recv_NPatchTotal != Recv_NPatchTotal_Check )
      Aux_Error( ERROR_INFO, "lv %d, Recv_NPatchTotal (%d) != expected (%d) !!\n",
                 lv, Recv_NPatchTotal, Recv_NPatchTotal_Check );
#  endif // #ifdef DEBUG_PARTICLE


// must NOT call "return" here even if Send/Recv_NPatchTotal==0 since this rank still needs to call Par_LB_SendParticleData()
// if ( Send_NPatchTotal == 0 )   return;
// if ( Recv_NPatchTotal == 0 )   return;


// 1. get the number of particles to be sent
   int *SendBuf_NParEachPatch = new int [Send_NPatchTotal];

   int PID, NParThisPatch, NSendParTotal = 0;

// loop over all target send patches
   for (int t=0; t<Send_NPatchTotal; t++)
   {
      PID           = Send_PIDList[t];
      NParThisPatch = amr->patch[0][lv][PID]->NPar;

#     ifdef DEBUG_PARTICLE
      if ( NParThisPatch < 0 )
         Aux_Error( ERROR_INFO, "NParThisPatch (%d) has not been calculated (lv %d, PID %d) !!\n",
                    NParThisPatch, lv, PID );
#     endif

      NSendParTotal            += NParThisPatch;
      SendBuf_NParEachPatch[t]  = NParThisPatch;
   } // for (int t=0; t<Send_NPatchTotal; t++)


// 2. prepare the particle data to be sent (and then remove these particles from this rank)
   const bool RemoveAllPar_Yes = true;

// reuse the MPI send buffer declared in LB_GetBufferData for better MPI performance
   real *SendBuf_ParDataEachPatch = LB_GetBufferData_MemAllocate_Send( NSendParTotal*PAR_NATT_TOTAL );

   real *SendPtr = SendBuf_ParDataEachPatch;
   long *ParList = NULL;
   long  ParID;

   for (int t=0; t<Send_NPatchTotal; t++)
   {
      PID           = Send_PIDList         [t];
      NParThisPatch = SendBuf_NParEachPatch[t];

//    skip patches with no particles
      if ( NParThisPatch == 0 )  continue;

      ParList = amr->patch[0][lv][PID]->ParList;

#     ifdef DEBUG_PARTICLE
      if ( ParList == NULL )
         Aux_Error( ERROR_INFO, "ParList == NULL for NParThisPatch (%d) > 0 (lv %d, PID %d) !!\n",
                    NParThisPatch, lv, PID );
#     endif

      for (int p=0; p<NParThisPatch; p++)
      {
         ParID = ParList[p];

//       2-1. store particle data into the MPI send buffer
         for (int v=0; v<PAR_NATT_TOTAL; v++)   *SendPtr++ = amr->Par->Attribute[v][ParID];

//       2-2. remove this particle from the particle repository of this rank
         amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_MPI );
      }

//    2-3. remove all particles in this send patch
      amr->patch[0][lv][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[lv], RemoveAllPar_Yes );
   } // for (int t=0; t<Send_NPatchTotal; t++)


// 3. send the number of particles and their attributes
   const bool Exchange_NPatchEachRank_No   = false;
   const bool Exchange_LBIdxEachRank_No    = false;
   const bool Exchange_ParDataEachRank_Yes = true;

   int  *SendBuf_NPatchEachRank   = Send_NPatchEachRank;
   int  *RecvBuf_NPatchEachRank   = Recv_NPatchEachRank;
   int  *RecvBuf_NParEachPatch    = NULL;    // will be allocated by Par_LB_SendParticleData and must be free'd later
   real *RecvBuf_ParDataEachPatch = NULL;    // a pointer to the MPI recv buffer declared in LB_GetBufferData
                                             // --> don't have to be free'd here

   long *SendBuf_LBIdxEachRank    = NULL;    // useless and does not need to be allocated
   long *RecvBuf_LBIdxEachRank    = NULL;    // useless and will not be allocated by Par_LB_SendParticleData

   int NRecvPatchTotal, NRecvParTotal;       // returned from Par_LB_SendParticleData

// note that we don't exchange NPatchEachRank (which is already known) and LBIdxEachRank (which is useless here)
   Par_LB_SendParticleData(
      PAR_NATT_TOTAL,
      SendBuf_NPatchEachRank, SendBuf_NParEachPatch, SendBuf_LBIdxEachRank, SendBuf_ParDataEachPatch, NSendParTotal,
      RecvBuf_NPatchEachRank, RecvBuf_NParEachPatch, RecvBuf_LBIdxEachRank, RecvBuf_ParDataEachPatch,
      NRecvPatchTotal, NRecvParTotal, Exchange_NPatchEachRank_No, Exchange_LBIdxEachRank_No, Exchange_ParDataEachRank_Yes,
      Timer, Timer_Comment );

#  ifdef DEBUG_PARTICLE
   if ( NRecvPatchTotal != Recv_NPatchTotal )
      Aux_Error( ERROR_INFO, "Total number of received patches (%d) != expected (%d) !!\n",
                 NRecvPatchTotal, Recv_NPatchTotal );
#  endif

// free the send buffer in advance to save memory
   delete [] SendBuf_NParEachPatch;


// 4. store the received particle data to the particle repository and link to each recv patch
   const real *RecvPtr = RecvBuf_ParDataEachPatch;

// 4-1. get the maximum number of particles in one recv patch
   int   NParThisPatch_Max;
   long *NewParIDList = NULL;

   NParThisPatch_Max = 0;
   for (int t=0; t<Recv_NPatchTotal; t++)    NParThisPatch_Max = MAX( NParThisPatch_Max, RecvBuf_NParEachPatch[t] );

   NewParIDList = new long [NParThisPatch_Max];


   for (int t=0; t<Recv_NPatchTotal; t++)
   {
      NParThisPatch = RecvBuf_NParEachPatch[t];

//    4-2. add particles to the particle repository
      for (int p=0; p<NParThisPatch; p++)
      {
         ParID    = amr->Par->AddOneParticle( RecvPtr );
         RecvPtr += PAR_NATT_TOTAL;

//       store the new particle index
         NewParIDList[p] = ParID;

//       we do not transfer inactive particles
#        ifdef DEBUG_PARTICLE
         if ( amr->Par->Attribute[PAR_MASS][ParID] < (real)0.0 )
            Aux_Error( ERROR_INFO, "Find inactive particle (ParID %d, Mass %14.7e) !!\n",
                       ParID, amr->Par->Attribute[PAR_MASS][ParID] );
#        endif
      }

//    4-3. add particles to the recv patch
      PID = Recv_PIDList[t];

      const real *PType = amr->Par->Type;
#     ifdef DEBUG_PARTICLE
//    do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//    may change after calling amr->Par->AddOneParticle
      const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      char Comment[100];
      sprintf( Comment, "%s", __FUNCTION__ );

      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParIDList, &amr->Par->NPar_Lv[lv],
                                           PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#     else
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParIDList, &amr->Par->NPar_Lv[lv], 
                                           PType );
#     endif
   } // for (int t=0; t<Recv_NPatchTotal; t++)


// 5. free memory
   delete [] RecvBuf_NParEachPatch;
   delete [] NewParIDList;

} // FUNCTION : Par_LB_ExchangeParticleBetweenPatch



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
