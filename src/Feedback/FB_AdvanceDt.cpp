#include "GAMER.h"

#ifdef FEEDBACK



// prototypes of built-in feedbacks
void FB_SNe( const int lv, const double TimeNew, const double TimeOld, const double dt,
             const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
             real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
             const int TID, RandomNumber_t *RNG );


// user-specified feedback to be set by a test problem initializer
void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                     const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                     real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
                     const int TID, RandomNumber_t *RNG ) = NULL;


// random number generators
extern RandomNumber_t *FB_RNG;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_AdvanceDt
// Description :  Feedback from particles to grids
//
// Note        :  1. Invoked by EvolveLevel()
//                2. FB_LEVEL must equal MAX_LEVEL for now
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                SaveSg_Flu : Sandglass to store the updated fluid data
//                SaveSg_Mag : Sandglass to store the updated B field
//
// Return      :  Update both grids and particles
//-------------------------------------------------------------------------------------------------------
void FB_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag )
{

// only work on FB_LEVEL for now
   if ( lv != FB_LEVEL )   return;


// nothing to do if all feedback options are disabled
   if ( ! FB_Any )   return;


// 1. collect particles for the sibling buffer patches
//    --> exclude father-sibling buffer patches (FaSibBufPatch_No) since currently we assume FB_LEVEL == MAX_LEVEL
//    --> for simplicity, we disable particle postition prediction (PredictPos_No) currently even though particles
//        just crossing from coarse (lv-1) to fine (lv) grids may have time greater than other particles at lv
//        --> to fix it, we will need to correct other particle attributes such as velocity too
   const bool TimingSendPar_Yes = true;
   const bool JustCountNPar_No  = false;
   const bool PredictPos_No     = false;
   const bool SibBufPatch_Yes   = true;
   const bool FaSibBufPatch_No  = false;
//###OPTIMIZATION: only collect necessary particle attributes
   const long ParAttBitIdx_In   = _PAR_TOTAL;

   Par_CollectParticle2OneLevel( lv, ParAttBitIdx_In, PredictPos_No, TimeNew, SibBufPatch_Yes, FaSibBufPatch_No,
                                 JustCountNPar_No, TimingSendPar_Yes );

#  ifdef DEBUG_PARTICLE
   if (  ! ( ParAttBitIdx_In & _PAR_POSX )  ||
         ! ( ParAttBitIdx_In & _PAR_POSY )  ||
         ! ( ParAttBitIdx_In & _PAR_POSZ )  )
      Aux_Error( ERROR_INFO, "ParAttBitIdx_In must include particle positions !!\n" );
#  endif



// 2. allocate a temporary particle repository to store the updated particle data
//    --> must initialize it since we will replace amr->Par->Attribute[] by ParAtt_Updated[]
//        at the end of this routine
//        --> otherwise, the data of particles skipped by this routine (e.g., those not on FB_LEVEL) will be incorrect
//        --> also to retain the information of inactive particles
//###OPTIMIZATION: only store the attributes being updated
//###OPTIMIZATION: only count particles on FB_LEVEL
   real *ParAtt_Updated[PAR_NATT_TOTAL];
   long  ParAttBitIdx_Out = _PAR_TOTAL;;

// do not update particle positions and accelerations
   ParAttBitIdx_Out &= ~( _PAR_POSX | _PAR_POSY | _PAR_POSZ );
#  ifdef STORE_PAR_ACC
   ParAttBitIdx_Out &= ~( _PAR_ACCX | _PAR_ACCY | _PAR_ACCZ );
#  endif

   for (int v=0; v<PAR_NATT_TOTAL; v++) {
      if ( ParAttBitIdx_Out & BIDX(v) ) {
         ParAtt_Updated[v] = new real [ amr->Par->ParListSize ];
         memcpy( ParAtt_Updated[v], amr->Par->Attribute[v], amr->Par->ParListSize*sizeof(real) );
      }
      else
         ParAtt_Updated[v] = NULL;
   }



// get the sibling index differences along different directions
   int NSibPID_Delta[26], *SibPID_Delta[26];

   TABLE_GetSibPID_Delta( NSibPID_Delta, SibPID_Delta );


// start of OpenMP parallel region
#  pragma omp parallel
   {

// thread ID
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif


// array to store the input and output fluid data
   real (*fluid_PG)[PS2][PS2][PS2] = new real [NCOMP_TOTAL][PS2][PS2][PS2];


// iterate over all real patches
#  pragma omp for schedule( runtime )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      const double PGCenter[3] = { amr->patch[0][lv][PID0+7]->EdgeL[0],
                                   amr->patch[0][lv][PID0+7]->EdgeL[1],
                                   amr->patch[0][lv][PID0+7]->EdgeL[2] };


//    3. prepare the fluid data to be updated
//    --> exclude magnetic field for now
//    --> use patch group as the basic unit
      const int  GhostZone_No        = 0;
      const int  NPG                 = 1;
#     ifndef MHD
      const int  OPT__MAG_INT_SCHEME = INT_NONE;
#     endif
      const bool IntPhase_No         = false;
      const real MinDens_No          = -1.0;
      const real MinPres_No          = -1.0;
      const real MinTemp_No          = -1.0;
      const bool DE_Consistency_No   = false;
//###OPTIMIZATION: only prepare the necessary fluid fields
      const long FluidBitIdx         = _TOTAL;

      Prepare_PatchData( lv, TimeNew, fluid_PG[0][0][0], NULL, GhostZone_No, NPG, &PID0, FluidBitIdx, _NONE,
                         OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_00, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );



//    4. collect patch information
      const int NNearbyPatchMax = 64;  // maximum number of neaby patches of a patch group (including 8 local patches)
      int Nearby_PID_List[NNearbyPatchMax], NNearbyPatch, SibPID0_List[26];

//    4-1. get nearby patches
      NNearbyPatch = 0;

//    local patches
      for (int PID=PID0; PID<PID0+8; PID++)  Nearby_PID_List[ NNearbyPatch ++ ] = PID;

//    sibling patches
      TABLE_GetSibPID_Based( lv, PID0, SibPID0_List );

//###OPTIMIZATION: skip sibling patches if the maximum feedback radius is zero
      for (int s=0; s<26; s++)
      {
         const int SibPID0 = SibPID0_List[s];   // first target patch in the sibling patch group

//       only consider leaf patches on FB_LEVEL (including both real and buffer patches)
         if ( SibPID0 >= 0 )
         for (int c=0; c<NSibPID_Delta[s]; c++)
         {
            const int SibPID = SibPID0 + SibPID_Delta[s][c];
            Nearby_PID_List[ NNearbyPatch ++ ] = SibPID;
         }
      }


//    4-2. record the coarse-fine boundaries
//         --> regard non-periodic boundaries as coarse-fine boundaries too
      bool CoarseFine[26];
      for (int s=0; s<26; s++)   CoarseFine[s] = ( SibPID0_List[s] < 0 ) ? true : false;



//    5. allocate arrays to store the local particle data
//       --> **local** means particles that affect the target patch group
//       --> include particles in both real and buffer patches
//       --> allocate the **maximum** required size among all nearby patches of a given patch group **just once**
//           for better performance
      real *ParAtt_Local[PAR_NATT_TOTAL];
      int  *ParSortID = NULL;
      int   NParMax   = -1;

      for (int t=0; t<NNearbyPatch; t++)
      {
         const int PID = Nearby_PID_List[t];

//       check both NPar and NPar_Copy (NPar_Copy may be -1, which is fine)
         NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar      );
         NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar_Copy );
      }

      if ( NParMax > 0 )
      {
         for (int v=0; v<PAR_NATT_TOTAL; v++)
         {
            if ( ParAttBitIdx_In & BIDX(v) )    ParAtt_Local[v] = new real [NParMax];
            else                                ParAtt_Local[v] = NULL;
         }

         ParSortID = new int [NParMax];   // it will fail if "long" is actually required for NParMax
      }


//    iterate over all nearby patches of the target patch group to apply feedback
      for (int t=0; t<NNearbyPatch; t++)
      {
         const int PID = Nearby_PID_List[t];


//       6. prepare the input particle data
//       6-1. get the particle list of the target patch
         long  *ParList = NULL;
         int    NPar;
         bool   UseParAttCopy;

//       the check "son == -1" is actually useless for now since we fix FB_LEVEL == MAX_LEVEL
         if ( amr->patch[0][lv][PID]->son == -1  &&  PID < amr->NPatchComma[lv][1] )
         {
            NPar          = amr->patch[0][lv][PID]->NPar;
            ParList       = amr->patch[0][lv][PID]->ParList;
            UseParAttCopy = false;

#           ifdef DEBUG_PARTICLE
            if ( amr->patch[0][lv][PID]->NPar_Copy != -1 )
               Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                          lv, PID, amr->patch[0][lv][PID]->NPar_Copy );
#           endif
         }

         else
         {
//          note that amr->patch[0][lv][PID]->NPar>0 is still possible
            NPar          = amr->patch[0][lv][PID]->NPar_Copy;
#           ifdef LOAD_BALANCE
            ParList       = NULL;
            UseParAttCopy = true;
#           else
            ParList       = amr->patch[0][lv][PID]->ParList_Copy;
            UseParAttCopy = false;
#           endif
         } // if ( amr->patch[0][lv][PID]->son == -1  &&  PID < amr->NPatchComma[lv][1] ) ... else ...


//       6-2. copy the particle data
//            --> we don't want to modify the input particle data during the iteration of different patch groups
//                since different patch groups will be affected by the same particles when feedback is non-local
//            --> if we modify the input particle data here, some patch groups may read the **updated** particle data
//            --> to solve this problem, we will store the updated particle data in a temporary particle repository
//                ParAtt_Updated[]
#        ifdef LOAD_BALANCE
         if ( UseParAttCopy ) {
            for (int v=0; v<PAR_NATT_TOTAL; v++) {
               if ( ParAttBitIdx_In & BIDX(v) ) {

#                 ifdef DEBUG_PARTICLE
                  if ( NPar > 0  &&  amr->patch[0][lv][PID]->ParAtt_Copy[v] == NULL )
                     Aux_Error( ERROR_INFO, "ParAtt_Copy == NULL for NPar (%d) > 0 (lv %d, PID %d, v %d) !!\n",
                                NPar, lv, PID, v );
#                 endif

                  for (int p=0; p<NPar; p++)
                     ParAtt_Local[v][p] = amr->patch[0][lv][PID]->ParAtt_Copy[v][p];
         }}}

         else
#        endif // #ifdef LOAD_BALANCE
         {
#           ifdef DEBUG_PARTICLE
            if ( NPar > 0  &&  ParList == NULL )
               Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PID %d) !!\n",
                          NPar, lv, PID );
#           endif

            for (int v=0; v<PAR_NATT_TOTAL; v++) {
               if ( ParAttBitIdx_In & BIDX(v) )
                  for (int p=0; p<NPar; p++)
                     ParAtt_Local[v][p] = amr->Par->Attribute[v][ ParList[p] ];
            }
         } // if ( UseParAttCopy ) ... else ...


//       6-3. sort particles by positions to fix their order
//            --> necessary when feedback involves random numbers
//            --> otherwise, the same particles accessed by different patches may have different random numbers
         Par_SortByPos( NPar, ParAtt_Local[PAR_POSX], ParAtt_Local[PAR_POSY], ParAtt_Local[PAR_POSZ], ParSortID );



//       6-4. periodic boundary conditions
//            --> assuming there are at least TWO patch groups along each spatial direction (i.e., NX0_TOT[*]/PS2>=2)
//            --> so each particle won't affect the same patch group more than once
         const bool   Periodic [3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                       OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                       OPT__BC_FLU[4] == BC_FLU_PERIODIC };
         const double HalfBox  [3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
         const int    ParPosIdx[3] = { PAR_POSX, PAR_POSY, PAR_POSZ };

         for (int d=0; d<3; d++) {
            if ( Periodic[d] ) {
               for (int p=0; p<NPar; p++) {

                  real *ParPos = ParAtt_Local[ ParPosIdx[d] ] + p;
                  const double dr = *ParPos - PGCenter[d];

                  if      ( dr > +HalfBox[d] )  *ParPos -= amr->BoxSize[d];
                  else if ( dr < -HalfBox[d] )  *ParPos += amr->BoxSize[d];
               }
         }} // for d, if ( Periodic[d] )



//       7. feedback
//       7-1. set the random seed
//            --> to get deterministic and different random numbers for all patches, reset the random seed of
//                each patch according to its location and counter
//            --> factors 1e2 and  1e8 are to make random seeds more different
         const long RSeed = FB_RSEED + amr->patch[0][lv][PID]->LB_Idx*100L + AdvanceCounter[lv]*100000000L;
         FB_RNG->SetSeed( TID, RSeed );


//       7-2. invoke all feedback routines
         if ( FB_SNE  )    FB_SNe     ( lv, TimeNew, TimeOld, dt, NPar, ParSortID, ParAtt_Local, fluid_PG,
                                        amr->patch[0][lv][PID0]->EdgeL, amr->dh[lv], CoarseFine, TID, FB_RNG );

         if ( FB_USER )    FB_User_Ptr( lv, TimeNew, TimeOld, dt, NPar, ParSortID, ParAtt_Local, fluid_PG,
                                        amr->patch[0][lv][PID0]->EdgeL, amr->dh[lv], CoarseFine, TID, FB_RNG );



//       8. store the updated particle data in ParAtt_Updated[]
//          --> only for particles in the central 8 patches
//              --> particles in the sibling patches will be updated and stored when applying feedback to these patches
//              --> avoid duplicate updates
//          --> different OpenMP threads work on different patch groups and thus won't update the same particles
         if ( t < 8  &&  amr->patch[0][lv][PID]->son == -1 ) {
            for (int v=0; v<PAR_NATT_TOTAL; v++)
               if ( ParAttBitIdx_Out & BIDX(v) )
                  for (int p=0; p<NPar; p++)
                     ParAtt_Updated[v][ ParList[p] ] = ParAtt_Local[v][p];
         }
      } // for (int t=0; t<NNearbyPatch; t++)



//    9. store the updated fluid data
//       --> different OpenMP threads work on different patch groups and thus won't update the same fluid data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID    = PID0 + LocalID;
         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 );
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );

         for (int v=0; v<NCOMP_TOTAL; v++)   {

            if ( FluidBitIdx & BIDX(v) )
            for (int k_o=0; k_o<PS1; k_o++)  {  const int k_i = Disp_k + k_o;
            for (int j_o=0; j_o<PS1; j_o++)  {  const int j_i = Disp_j + j_o;
            for (int i_o=0; i_o<PS1; i_o++)  {  const int i_i = Disp_i + i_o;

               amr->patch[SaveSg_Flu][lv][PID]->fluid[v][k_o][j_o][i_o] = fluid_PG[v][k_i][j_i][i_i];

            }}}
         }
      } // for (int LocalID=0; LocalID<8; LocalID++)


//    free memory
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParAtt_Local[v];
      delete [] ParSortID;

   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

// free memory
   delete [] fluid_PG;

   } // end of OpenMP parallel region



// 10. store the updated particle data
   for (int v=0; v<PAR_NATT_TOTAL; v++) {
      if ( ParAttBitIdx_Out & BIDX(v) )
         memcpy( amr->Par->Attribute[v], ParAtt_Updated[v], amr->Par->ParListSize*sizeof(real) );
   }



// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParAtt_Updated[v];
   Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_Yes, FaSibBufPatch_No );

} // FUNCTION : FB_AdvanceDt



#endif // #ifdef FEEDBACK
