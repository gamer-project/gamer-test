#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_CorrAfterAllSync
// Description :  Apply various corrections after all levels are synchronized
//
// Note        :  1. Corrections include:
//                   1-1. Synchronize all particles
//                        --> Since particles just crossing coarse-fine boundaries may not be synchronized
//                            with the grid data
//                        --> Must work with STORE_PAR_ACC
//                   1-2. Restrict data
//                        --> This is mainly for bitwise reproducibility when restarting simulations from
//                            C binary outputs
//                            --> Because C binary outputs do not store data in non-leaf patches, we must
//                                apply the restrict operation to obtain these data during restart. But for
//                                father patches with new sons, the round-off erros for these father patches
//                                in the original simulations and in the restart process can be different.
//                   1-3. Recalculate gravitational potential
//                        --> This is for both improving accuracy and bitwise reproducibility during restart
//                            --> Accuracy: because it uses newly restricted density to recalculate potential
//                            --> Reproducibility: because we always recalculate potential during restart
//                   1-4. Recalculate particle acceleration
//                        --> Because particles have been synchronized and potential has been recalculated
//                        --> Must work with STORE_PAR_ACC
//                2. Invoked by main() and Output_DumpData()
//                3. Do nothing if there are only base-level patches
//                4. Strictly speaking, this function belongs to both Fluid, SelfGravity, and Particle categories
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Flu_CorrAfterAllSync()
{

// nothing to do if there are only base-level patches
   if ( NLEVEL == 1  ||  NPatchTotal[1] == 0 )     return;


   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "   %s                     ...\n", __FUNCTION__ );


// 1. synchronize all particles
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC )
   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "      synchronize particles                 ... " );

   if (  Par_Synchronize( Time[0], PAR_SYNC_FORCE ) != 0  )
      Aux_Error( ERROR_INFO, "particle synchronization failed !!\n" );

// particles may cross patch boundaries after synchronization
   const bool TimingSendPar_No = false;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Par_PassParticle2Sibling( lv, TimingSendPar_No );
      Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_No, NULL_INT, NULL );
   }

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#  endif


// 2. restrict data
   if ( OPT__FIXUP_RESTRICT )
   for (int lv=MAX_LEVEL-1; lv>=0; lv--)
   {
      if ( NPatchTotal[lv+1] == 0 )    continue;

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "      restrict data at Lv %2d                ... ", lv );

//    we do not restrict potential since it will be recalculated anyway
      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT, _TOTAL, _MAG );

#     ifdef LOAD_BALANCE
      LB_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT, _TOTAL, _MAG, NULL_INT );
#     endif

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_AFTER_FIXUP, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }


// 3. recalculate gravitational potential
#  ifdef GRAVITY
   if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
   for (int lv=0; lv<=MAX_LEVEL; lv++)
   {
      if ( NPatchTotal[lv] == 0 )   break;

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "      recalculate potential at Lv %2d        ... ", lv );

      if ( lv > 0 )
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Rho_ParaBuf, USELB_YES );

      Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false, false );

      if ( lv > 0 )
      Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

   } // for (int lv=0; lv<=MAX_LEVEL; lv++) if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
#  endif // #ifdef GRAVITY


// 4. recalculate particle acceleration
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC  && defined GRAVITY )
   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "      recalculate particle acceleration     ... " );

   const bool StoreAcc_Yes    = true;
   const bool UseStoredAcc_No = false;

   for (int lv=0; lv<NLEVEL; lv++)
   Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#  endif


   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "   %s                     ... done\n", __FUNCTION__ );

} // FUNCTION : Flu_CorrAfterAllSync
