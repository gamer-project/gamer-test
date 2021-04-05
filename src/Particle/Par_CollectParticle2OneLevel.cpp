#include "GAMER.h"

#ifdef PARTICLE

# ifndef LOAD_BALANCE
static void CollectParticle( const int FaLv, const int FaPID, int &NPar_SoFar, long *ParList );
#endif

// flag (declared in Prepare_PatchData.cpp) for checking whether Par_CollectParticle2OneLevel() has been called before
// preparing either _PAR_DENS or _TOTAL_DENS data in Prepare_PatchData()
extern bool Particle_Collected;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CollectParticle2OneLevel
// Description :  Collect particles from all descendants (sons, grandsons, ...) to patches at a target level
//                --> Will call Par_LB_CollectParticle2OneLevel for LOAD_BALANCE
//
// Note        :  1. ParList in all descendants will NOT be changeed after calling this function
//                2. This function assumes that father and all descendants are always in the same rank
//                   --> Does NOT work with LOAD_BALANCE, and it will call Par_LB_CollectParticle2OneLevel for that
//                3. The array "ParList_Copy" will be allocated for the target patch
//                   --> Must be deleted afterward by calling Par_CollectParticle2OneLevel_FreeMemory
//                       (e.g., in the function "Gra_AdvanceDt")
//                4. Only work on non-leaf patches
//                   --> For leaf patches the ParList_Copy will NOT be allocated
//                   --> However, it does take into account the occasional situation where non-leaf patches may
//                       have NPar>0 temporarily after updating particle position.
//                       It's because particles travelling from coarse to fine grids will stay in coarse grids
//                       temporarily until the velocity correction is done.
//                       --> For these patches, NPar_Copy will be **the sum of NPar and the number of particles
//                           collected from other patches**, and ParList_Copy (or ParAtt_Copy) will contain
//                           information of particles belonging to NPar as well.
//                       --> It makes implementation simplier. **For leaf real patches, one only needs to consider
//                           NPar and ParList. While for all other patches, one only needs to consider NPar_Copy,
//                           ParList_Copy (or ParAtt_Copy). One never needs to consider both.**
//                5. When using OpenMP, one must ensure that different threads do NOT invoke this function
//                   for the same patch at the same time !!!
//                   --> Because this function will modify "NPar_Copy & ParList_Copy" for the target patch
//                6. Invoked by Gra_AdvanceDt(), Flag_Real(), Output_DumpData_Total(), Output_DumpData_Total_HDF(),
//                   and Flu_CorrAfterAllSync()
//                7. When turning on SibBufPatch in LOAD_BALANCE, this function (which will call
//                   Par_LB_CollectParticle2OneLevel) will also collect particles for sibling-buffer patchesat FaLv
//                   --> Moreover, if FaSibBufPatch is also on, it will also collect particles for
//                       father-sibling-buffer patches at FaLv-1 (if FaLv > 0)
//                       --> Useful for constructing the density field at FaLv for the Poisson solver at FaLv
//                8. Option "JustCountNPar" can be used to count the number of particles in each real patch at FaLv
//                   --> Do NOT collect particle indices
//                       --> ParList_Copy will NOT be allocated
//                   --> Particle count is stored in NPar_Copy
//
// Parameter   :  FaLv          : Target refinement leve
//                AttBitIdx     : Bitwise indices of the target particle attributes (e.g., _PAR_MASS | _PAR_VELX)
//                                --> A user-defined attribute with an integer index AttIntIdx returned by
//                                    AddParticleAttribute() can be converted to a bitwise index by BIDX(AttIntIdx)
//                                --> Used by LOAD_BALANCE only
//                PredictPos    : true --> Predict particle position to TargetTime (for LOAD_BALANCE only)
//                TargetTime    : Target time for predicting the particle position (for LOAD_BALANCE only)
//                SibBufPatch   : true --> Collect particles for sibling-buffer patches at FaLv as well
//                                         (for LOAD_BALANCE only)
//                FaSibBufPatch : true --> Collect particles for father-sibling-buffer patches at FaLv-1 as well
//                                         (do nothing if FaLv==0) (for LOAD_BALANCE only)
//                JustCountNPar : Just count the number of particles in each real patch at FaLv. Don't collect
//                                particle indices (or collect particle mass and position for LOAD_BALANCE)
//                TimingSendPar : Measure the elapsed time of Par_LB_SendParticleData() in
//                                Par_LB_CollectParticle2OneLevel() (for LOAD_BALANCE only)
//
// Return      :  NPar_Copy and ParList_Copy (if JustCountNPar == false) for all non-leaf real patches at FaLv
//-------------------------------------------------------------------------------------------------------
void Par_CollectParticle2OneLevel( const int FaLv, const long AttBitIdx, const bool PredictPos, const double TargetTime,
                                   const bool SibBufPatch, const bool FaSibBufPatch, const bool JustCountNPar,
                                   const bool TimingSendPar )
{

// set this flag to true to indicate that this function has been called
// --> must be set before invoking the load-balance alternative routine Par_LB_CollectParticle2OneLevel()
   Particle_Collected = true;


// call the parallel version instead
#  ifdef LOAD_BALANCE
// note that if SibBufPatch or FaSibBufPatch is on, we need to call Par_LB_CollectParticle2OneLevel()
// even when FaLv == MAX_LEVEL
   Par_LB_CollectParticle2OneLevel( FaLv, AttBitIdx, PredictPos, TargetTime, SibBufPatch, FaSibBufPatch, JustCountNPar, TimingSendPar );

   return;


#  else


// nothing to do for the max level
   if ( FaLv >= MAX_LEVEL )   return;


//###NOTE: OpenMP may not improve performance here
#  pragma omp parallel for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
//    check
#     ifdef DEBUG_PARTICLE
      if ( amr->patch[0][FaLv][FaPID]->NPar_Copy != -1  ||  amr->patch[0][FaLv][FaPID]->ParList_Copy != NULL )
         Aux_Error( ERROR_INFO, "Some variables have been initialized already (FaLv %d, FaPID %d, NPar_Copy %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Copy );
#     endif


//    nothing to do if this father has no son
//    --> leaf patches will always have NPar_Copy == -1
      if ( amr->patch[0][FaLv][FaPID]->son == -1 )    continue;


//    1. get the total number of particles in all descendants
      amr->patch[0][FaLv][FaPID]->NPar_Copy = Par_CountParticleInDescendant( FaLv, FaPID );


//    2. add the number of particles temporarily residing in this patch waiting for the velocity correction in KDK
      amr->patch[0][FaLv][FaPID]->NPar_Copy += amr->patch[0][FaLv][FaPID]->NPar;

#     ifdef DEBUG_PARTICLE
//    check if these particles are indeed waiting for the velocity correction (i.e., ParTime = -dt_half < 0.0 for KDK)
      if ( amr->Par->Integ == PAR_INTEG_KDK )
      for (int p=0; p<amr->patch[0][FaLv][FaPID]->NPar; p++)
      {
         const long ParID = amr->patch[0][FaLv][FaPID]->ParList[p];

         if ( amr->Par->Time[ParID] >= (real)0.0 )
            Aux_Error( ERROR_INFO, "This particle shouldn't be here (FaLv %d, FaPID %d, ParID %ld, ParTime %21.14e) !!\n",
                       FaLv, FaPID, ParID, amr->Par->Time[ParID] );
      }
#     endif


//    3. nothing to do if we only want to count particles (or if this patch has no particles at all)
      if ( JustCountNPar  ||  amr->patch[0][FaLv][FaPID]->NPar_Copy == 0 )    continue;


//    4. allocate the array ParList_Copy
      amr->patch[0][FaLv][FaPID]->ParList_Copy = new long [ amr->patch[0][FaLv][FaPID]->NPar_Copy ];


//    5. gather the particle lists from all descendants
      int NPar_SoFar = 0;
      CollectParticle( FaLv, FaPID, NPar_SoFar, amr->patch[0][FaLv][FaPID]->ParList_Copy );


//    6. add particles temporarily residing in this patch
      for (int p=0; p<amr->patch[0][FaLv][FaPID]->NPar; p++)
      {
         const long ParID = amr->patch[0][FaLv][FaPID]->ParList[p];

         amr->patch[0][FaLv][FaPID]->ParList_Copy[ NPar_SoFar ++ ] = ParID;
      }


//    7. check the particle count
#     ifdef DEBUG_PARTICLE
      if ( NPar_SoFar != amr->patch[0][FaLv][FaPID]->NPar_Copy )
         Aux_Error( ERROR_INFO, "NPar_SoFar (%d) != NPar_Copy (%d) !!\n",
                    NPar_SoFar, amr->patch[0][FaLv][FaPID]->NPar_Copy );
#     endif
   } // for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)


#  endif // #ifdef LOAD_BALANCE ... else ...


} // FUNCTION : Par_CollectParticle2OneLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CollectParticle2OneLevel_FreeMemory
// Description :  Release the memory allocated by Par_CollectParticle2OneLevel
//
// Note        :  1. Invoded by Gra_AdvanceDt (and Main when DEBUG is on)
//                2. For LOAD_BALANCE, this function will call the alternative function
//                   "Par_LB_CollectParticle2OneLevel_FreeMemory"
//
// Parameter   :  FaLv          : Target refinement leve
//                SibBufPatch   : true --> Release memory for sibling-buffer patches at FaLv as well (for LOAD_BALANCE only)
//                FaSibBufPatch : true --> Release memory for father-sibling-buffer patches at FaLv-1 as well
//                                         (do nothing if FaLv==0) (for LOAD_BALANCE only)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_CollectParticle2OneLevel_FreeMemory( const int FaLv, const bool SibBufPatch, const bool FaSibBufPatch )
{

// set this flag to false to indicate that Par_CollectParticle2OneLevel has NOT been called
// --> must be set before invoking the load-balance alternative routine "Par_LB_CollectParticle2OneLevel_FreeMemory"
   Particle_Collected = false;


#  ifdef LOAD_BALANCE

   Par_LB_CollectParticle2OneLevel_FreeMemory( FaLv, SibBufPatch, FaSibBufPatch );

#  else

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->ParList_Copy != NULL )
      {
         delete [] amr->patch[0][FaLv][FaPID]->ParList_Copy;
         amr->patch[0][FaLv][FaPID]->ParList_Copy = NULL;
      }

//    all patches should have NPar_Copy == -1 to indicate that it has not been calculated
      amr->patch[0][FaLv][FaPID]->NPar_Copy = -1;
   }

#  endif // #ifdef LOAD_BALANCE ... else ...

} // FUNCTION : Par_CollectParticle2OneLevel_FreeMemory



# ifndef LOAD_BALANCE
//-------------------------------------------------------------------------------------------------------
// Function    :  CollectParticle
// Description :  Collect particles from all descendants (sons, grandsons, ...) of the target patch
//
// Note        :  This function will search over all descendants recursively
//
// Parameter   :  FaLv       : Father patch level
//                FaPID      : Father patch ID
//                NPar_SoFar : Number of particles counted so for
//                             --> used as the starting array index in ParList
//                ParList    : Array to store the particle IDs
//
// Return      :  NPar_SoFar, ParList
//-------------------------------------------------------------------------------------------------------
void CollectParticle( const int FaLv, const int FaPID, int &NPar_SoFar, long *ParList )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );
#  endif


   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;

   if ( SonPID0 == -1 )    return;  // nothing to do if the target patch has no son
   else
   {
      for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      {
//       search over all grandsons recursively
         if ( amr->patch[0][SonLv][SonPID]->son != -1 )  CollectParticle( SonLv, SonPID, NPar_SoFar, ParList );

         else
            for (int p=0; p<amr->patch[0][SonLv][SonPID]->NPar; p++)
               ParList[ NPar_SoFar ++ ] = amr->patch[0][SonLv][SonPID]->ParList[p];
      }
   }

} // FUNCTION : CollectParticle
#endif // #ifndef LOAD_BALANCE



#endif // #ifdef PARTICLE
