#include "GAMER.h"

#ifdef PARTICLE

#ifdef TIMING
extern Timer_t *Timer_Par_MPI[NLEVEL][6];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Sibling
// Description :  Pass particles lying outside the current home patches to the sibling patches
//
// Note        :  1. We have assumed that all particles will only go to the ***NEARBY 26 patches***
//                   --> In the cases that particles cross the coarse-fine boundary, particles are sent to the
//                       corresponding sibling-son (cross C->F) or father-sibling (cross F->C) patches
//                   --> But note that for C->F we only send particles to the sibling patches and do NOT
//                       send them to the sibling-son patches in this function. These particles are stored in
//                       the sibling patches temporarily for the velocity correction in KDK (because we don't
//                       have potential at lv+1 level at this point) and will be sent to lv+1 level after that.
//                2. After calculating the target sibling patch, the particle position will be re-mapped
//                   into the simulation domain if periodic B.C. is assumed
//                3. Particles transferred to buffer patches (at either lv or lv-1) will be resent to their
//                   corresponding real patches by calling Par_LB_ExchangeParticleBetweenPatch()
//
// Parameter   :  lv            : Target refinement level
//                TimingSendPar : Measure the elapsed time of Par_LB_SendParticleData(), which is called by
//                                Par_LB_ExchangeParticleBetweenPatch()
//                                --> LOAD_BALANCE only
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Sibling( const int lv, const bool TimingSendPar )
{

   const int    FaLv             = lv - 1;
   const bool   RemoveAllPar_No  = false;
   const int    MirSib[26]       = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const double GuessEscapeRatio = 0.5*( CUBE(PS1) - CUBE(PS1-2) )/CUBE(PS1);
   const int    SibID[3][3][3]   = {  { {18, 10, 19}, {14,  4, 16}, {20, 11, 21} },
                                      { { 6,  2,  7}, { 0, -1,  1}, { 8,  3,  9} },
                                      { {22, 12, 23}, {15,  5, 17}, {24, 13, 25} }  };
   const double dh_min           = amr->dh[TOP_LEVEL];
   const double BoxEdge[3]       = { (NX0_TOT[0]*(1<<TOP_LEVEL))*dh_min,
                                     (NX0_TOT[1]*(1<<TOP_LEVEL))*dh_min,
                                     (NX0_TOT[2]*(1<<TOP_LEVEL))*dh_min }; // prevent from the round-off error problem
// ParPos should NOT be used after calling Par_LB_ExchangeParticleBetweenPatch() since amr->Par->Attribute may be reallocated
   real *ParPos[3]               = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *PType = amr->Par->Type;

   int     NPar_Remove_Tot=0;
   int     NPar, NGuess, NPar_Remove, ArraySize[26], ijk[3], Side, TSib, SibPID, FaPID, FaSib, FaSibPID;
   long    ParID;
   int    *RemoveParList;
   double *EdgeL, *EdgeR;

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

#  pragma omp parallel private( NPar, NGuess, NPar_Remove, ArraySize, ijk, TSib, ParID, RemoveParList, EdgeL, EdgeR )
   {

#  pragma omp for reduction( +:NPar_Remove_Tot ) schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
// loop over all **real** patches
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      NPar  = amr->patch[0][lv][PID]->NPar;
      EdgeL = amr->patch[0][lv][PID]->EdgeL;
      EdgeR = amr->patch[0][lv][PID]->EdgeR;


//    0. skip patches with no particles
      if ( NPar == 0 )  continue;


//    1. allocate memory for the escaping particle list
      NGuess = MAX( 1, int(NPar*GuessEscapeRatio/26) );

      for (int s=0; s<26; s++)
      {
         ArraySize                           [s] = NGuess;
         amr->patch[0][lv][PID]->NPar_Escp   [s] = 0;
         amr->patch[0][lv][PID]->ParList_Escp[s] = NULL;    // allocate only when necessary (==> better performance)
      }


//    2. record the escaping particles
      RemoveParList = new int [NPar];
      NPar_Remove   = 0;

      for (int p=0; p<NPar; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

         for (int d=0; d<3; d++)
         {
//          2-1. check if particles lie outside the patch
            ijk[d] = ( ParPos[d][ParID] < EdgeL[d] ) ? 0 : (ParPos[d][ParID] < EdgeR[d]) ? 1 : 2;

//          2-2. reset particle position for periodic B.C.
//               --> note that EdgeL/R in amr->patch always assumes periodicity
//               --> OK when calling patch->AddParticle() in the debug mode
            if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC  &&  ijk[d] != 1 )
            {
               if ( ParPos[d][ParID] < 0.0 )
               {
                  ParPos[d][ParID] = ParPos[d][ParID] + BoxEdge[d];

//                prevent from the round-off error problem
                  if ( ParPos[d][ParID] == BoxEdge[d] )
                  {
                     ijk[d]           = 1;
                     ParPos[d][ParID] = 0.0;
                  }

#                 ifdef DEBUG_PARTICLE
                  if ( ParPos[d][ParID] < 0.0  ||  ParPos[d][ParID] >= BoxEdge[d] )
                     Aux_Error( ERROR_INFO, "Dim %d: ParPos[%ld] (%20.14e) < 0.0  ||  >= BoxEdge (%20.14e) !!\n",
                                d, ParID, ParPos[d][ParID], BoxEdge[d] );
#                 endif
               }

               else if ( ijk[d] == 2  &&  ParPos[d][ParID] >= BoxEdge[d] )
               {
                  ParPos[d][ParID] = ParPos[d][ParID] - BoxEdge[d];

#                 ifdef DEBUG_PARTICLE
                  if ( BoxEdge[d] != EdgeR[d] )
                     Aux_Error( ERROR_INFO, "Dim %d: BoxEdge (%20.14e) != EdgeR (%20.14e) !!\n",
                                d, BoxEdge[d], EdgeR[d] );

                  if ( ParPos[d][ParID] < 0.0  ||  ParPos[d][ParID] >= BoxEdge[d] )
                     Aux_Error( ERROR_INFO, "Dim %d: ParPos[%ld] (%20.14e) < 0.0  ||  >= BoxEdge (%20.14e) !!\n",
                                d, ParID, ParPos[d][ParID], BoxEdge[d] );
#                 endif
               }
            } // if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC  &&  ijk[d] != 1 )
         } // for (int d=0; d<3; d++)

         TSib = SibID[ ijk[2] ][ ijk[1] ][ ijk[0] ];


//       2-3. remove particles lying outside the active region for non-periodic B.C. (by setting mass as PAR_INACTIVE_OUTSIDE)
         if (  !PeriodicAllDir  &&  !Par_WithinActiveRegion( ParPos[0][ParID], ParPos[1][ParID], ParPos[2][ParID] )  )
         {
            RemoveParList[ NPar_Remove ++ ] = p;
            NPar_Remove_Tot ++;

//          use OpenMP critical construct since RemoveOneParticle will modify NPar_Active/Inactive, which are global variables
//          --> note that the order of which thread calls RemoveOneParticle() is nondeterministic and may change from run to run
//              --> order of particles stored in the particle repository (i.e., their particle ID) may change from run to run
//              --> particle text file may change from run to run since it's dumped according to the order of particle ID
//          --> but it's not an issue since the actual data of each particle will not be affected
#           pragma omp critical
            amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_OUTSIDE );

            if ( OPT__VERBOSE )
               Aux_Message( stderr, "\nWARNING : removing particle %10d (Pos = [%14.7e, %14.7e, %14.7e], Time = %13.7e)\n",
                            ParID, ParPos[0][ParID], ParPos[1][ParID], ParPos[2][ParID], Time[lv] );
         }


//       2-4. deal with escaping particles (i.e., particles lying outside the patch but still within the active region)
         else if ( TSib != -1 )
         {
            RemoveParList[ NPar_Remove ++ ] = p;
            NPar_Remove_Tot ++;

//          allocate ParList_Escp if it's not allocated yet
            if ( amr->patch[0][lv][PID]->ParList_Escp[TSib] == NULL )
               amr->patch[0][lv][PID]->ParList_Escp[TSib] = (long*)malloc( ArraySize[TSib]*sizeof(long) );

//          make sure that there is enough memory
            if ( amr->patch[0][lv][PID]->NPar_Escp[TSib] >= ArraySize[TSib] )
            {
               ArraySize[TSib] += NGuess;
               amr->patch[0][lv][PID]->ParList_Escp[TSib] = (long*)realloc( amr->patch[0][lv][PID]->ParList_Escp[TSib],
                                                                            ArraySize[TSib]*sizeof(long) );
            }

            amr->patch[0][lv][PID]->ParList_Escp[TSib][ amr->patch[0][lv][PID]->NPar_Escp[TSib] ++ ] = ParID;

#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Mass[ParID] < 0.0 )
               Aux_Error( ERROR_INFO, "Storing escaping particles which have been removed (lv %d, PID %d, TSib %d, ParID %d) !!\n",
                          lv, PID, TSib, ParID );

//          usually it happens when particle position is NaN
            if ( amr->patch[0][lv][PID]->sibling[TSib] < -1 )
            {
               Aux_Message( stderr, "ERROR : This particle lies outside the simulation box (lv %d, PID %d, TSib %d, sib %d, ParID %ld) !!\n",
                            lv, PID, TSib, amr->patch[0][lv][PID]->sibling[TSib], ParID );
               Aux_Message( stderr, "        --> ParPos = (%21.14e, %21.14e, %21.14e)\n",
                            ParPos[0][ParID], ParPos[1][ParID], ParPos[2][ParID] );
#              ifdef MHD
               const int MagSg = amr->MagSg[lv];
#              else
               const int MagSg = NULL_INT;
#              endif
               Output_Patch( lv, PID, amr->FluSg[lv], amr->PotSg[lv], MagSg, "debug" );
               MPI_Exit();
            }
#           endif
         } // else if ( TSib != -1 )
      } // for (int p=0; p<NPar; p++)


//    3. remove the escaping particles (set amr->Par->NPar_Lv later due to OpenMP)
      amr->patch[0][lv][PID]->RemoveParticle( NPar_Remove, RemoveParList, NULL, RemoveAllPar_No );
      delete [] RemoveParList;

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

   } // end of OpenMP parallel region


// update the global variable NPar_Lv after the OpenMP parallel region
   amr->Par->NPar_Lv[lv] -= NPar_Remove_Tot;


   if ( NPar_Remove_Tot > 0 )
   {
//    loop over all real **and buffer** patches
      for (int PID=0; PID<amr->num[lv]; PID++)
      {
//       4. gather the escaping particles from the 26 sibling patches (coarse --> coarse)
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

//          SibPID can be negative for non-periodic BC.
            if ( SibPID >= 0  &&  amr->patch[0][lv][SibPID]->NPar_Escp[ MirSib[s] ] > 0 )
            {
//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
#              ifdef DEBUG_PARTICLE
               if ( SibPID >= amr->NPatchComma[lv][1] )
                  Aux_Error( ERROR_INFO, "buffer patch cannot have escaping particles (PID %d, s %d, SibPID %d, NPar_Escp %d) !!\n",
                             PID, s, SibPID, amr->patch[0][lv][SibPID]->NPar_Escp[ MirSib[s] ] );

               char Comment[100];
               sprintf( Comment, "%s C->C", __FUNCTION__ );
               amr->patch[0][lv][PID]->AddParticle( amr->patch[0][lv][SibPID]->NPar_Escp   [ MirSib[s] ],
                                                    amr->patch[0][lv][SibPID]->ParList_Escp[ MirSib[s] ],
                                                   &amr->Par->NPar_Lv[lv], PType,
                                                    (const real **)ParPos, amr->Par->NPar_AcPlusInac, Comment );
#              else
               amr->patch[0][lv][PID]->AddParticle( amr->patch[0][lv][SibPID]->NPar_Escp   [ MirSib[s] ],
                                                    amr->patch[0][lv][SibPID]->ParList_Escp[ MirSib[s] ],
                                                   &amr->Par->NPar_Lv[lv], PType );
#              endif
            }
         } // for (int s=0; s<26; s++)


//       5. for patches with sons, pass particles to their sons (coarse --> fine)
//       *** we now do this after the correction step of KDK so that particles just travel from lv to lv+1
//       *** can have their velocity corrected at lv first (because we don't have potential at lv+1 at this point)
//       if ( amr->patch[0][lv][PID]->son != -1 )  Par_PassParticle2Son_SinglePatch( lv, PID );
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


//###NOTE : NO OpenMP since particles from different patches can enter the same father-sibling patch
//    6. pass particles to the father-sibling patches (fine --> coarse)
      if ( lv > 0 )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int s=0; s<26; s++)
      {
         SibPID = amr->patch[0][lv][PID]->sibling[s];

//       do not pass particles if SibPID < -1 (--> outside the simulation box)
//       (actually we already guarantee that NPar_Escp = 0 if SibPID < -1)
         if ( SibPID == -1  &&  amr->patch[0][lv][PID]->NPar_Escp[s] > 0 )
         {
//          find the correct father->sibling patch index
            for (int d=0; d<3; d++)
            {
               ijk[d] = TABLE_01( s,     'x'+d, 0, 1, 2 );
               Side   = TABLE_02( PID%8, 'x'+d, -1, +1 );

               if (  ( Side < 0 && ijk[d] == 2 )  ||  ( Side > 0 && ijk[d] == 0 )  )   ijk[d] = 1;
            }

            FaSib = SibID[ ijk[2] ][ ijk[1] ][ ijk[0] ];
            FaPID = amr->patch[0][lv][PID]->father;

#           ifdef DEBUG_PARTICLE
            if ( FaSib == -1 )
               Aux_Error( ERROR_INFO, "FaSib == -1 (lv %d, PID %d, Sib %d) !!\n", lv, PID, s );

            if ( FaPID < 0 )
               Aux_Error( ERROR_INFO, "FaPID = %d < 0 (lv %d, PID %d) !!\n", FaPID, lv, PID );
#           endif

            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[FaSib];

#           ifdef DEBUG_PARTICLE
            if ( FaSibPID < 0 )
               Aux_Error( ERROR_INFO, "FaSibPID < 0 (FaLv %d, FaPID %d, Sib %d) !!\n", FaLv, FaPID, s );
#           endif


//          add particles to the target father->sibling patch (which can be real of buffer patches)
//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
#           ifdef DEBUG_PARTICLE
            char Comment[100];
            sprintf( Comment, "%s F->C", __FUNCTION__ );
            amr->patch[0][FaLv][FaSibPID]->AddParticle( amr->patch[0][lv][PID]->NPar_Escp   [s],
                                                        amr->patch[0][lv][PID]->ParList_Escp[s],
                                                       &amr->Par->NPar_Lv[FaLv], PType,
                                                        (const real **)ParPos, amr->Par->NPar_AcPlusInac, Comment );
#           else
            amr->patch[0][FaLv][FaSibPID]->AddParticle( amr->patch[0][lv][PID]->NPar_Escp   [s],
                                                        amr->patch[0][lv][PID]->ParList_Escp[s],
                                                       &amr->Par->NPar_Lv[FaLv], PType );
#           endif
         }
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++); for (int s=0; s<26; s++)
   } // if ( NPar_Remove_Tot > 0 )


// 7. send particles from buffer patches to the corresponding real patches
//    --> note that after calling the following rourtines, some particles may reside in **non-leaf** real patches
//    --> they will be sent again to leaf real patches after the velocity correction operation
//        --> by Par_PassParticle2Son_MultiPatch()
#  ifdef LOAD_BALANCE

   Timer_t *Timer[2] = { NULL, NULL };
   char Timer_Comment[2][20];

#  ifdef TIMING
   if ( TimingSendPar )
   {
      Timer[0] = Timer_Par_MPI[lv][0];
      Timer[1] = Timer_Par_MPI[lv][1];

      sprintf( Timer_Comment[0], "%3d %15s", lv, "Par_2Sib-Sib" );
      sprintf( Timer_Comment[1], "%3d %15s", lv, "Par_2Sib-FaSib" );
   }
#  endif

// 7-1. sibling-buffer patches at lv
   Par_LB_ExchangeParticleBetweenPatch(
      lv,
      amr->Par->B2R_Buff_NPatchTotal[lv][0], amr->Par->B2R_Buff_PIDList[lv][0], amr->Par->B2R_Buff_NPatchEachRank[lv][0],
      amr->Par->B2R_Real_NPatchTotal[lv][0], amr->Par->B2R_Real_PIDList[lv][0], amr->Par->B2R_Real_NPatchEachRank[lv][0],
      Timer[0], Timer_Comment[0] );

// 7-2. father-sibling-buffer patches at lv-1 (FaLv)
//      --> note that XXX[lv][1] is for exchanging patches at lv-1
   if ( FaLv >= 0 )
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      amr->Par->B2R_Buff_NPatchTotal[lv][1], amr->Par->B2R_Buff_PIDList[lv][1], amr->Par->B2R_Buff_NPatchEachRank[lv][1],
      amr->Par->B2R_Real_NPatchTotal[lv][1], amr->Par->B2R_Real_PIDList[lv][1], amr->Par->B2R_Real_NPatchEachRank[lv][1],
      Timer[1], Timer_Comment[1] );

// 7-3. check: no buffer patches at lv and lv-1 can have particles at this point
#  ifdef DEBUG_PARTICLE
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][3]; PID++)
      if ( amr->patch[0][lv][PID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "Buffer patch has particles (lv %d, PID %d, NPar %d) !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar );
   if ( FaLv >= 0 )
   for (int FaPID=amr->NPatchComma[FaLv][1]; FaPID<amr->NPatchComma[FaLv][3]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "Buffer patch has particles (FaLv %d, FaPID %d, NPar %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );
#  endif
#  endif // #ifdef LOAD_BALANCE


// 8. get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );


// 9. free memory
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int s=0; s<26; s++)
   {
      if ( amr->patch[0][lv][PID]->ParList_Escp[s] != NULL )   free( amr->patch[0][lv][PID]->ParList_Escp[s] );

      amr->patch[0][lv][PID]->ParList_Escp[s] = NULL;
      amr->patch[0][lv][PID]->NPar_Escp   [s] = -1;      // -1: indicate that it has not been calculated yet
   }

} // FUNCTION : Par_PassParticle2Sibling



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_WithinActiveRegion
// Description :  Check whether the input coordinates are within the active region
//
// Note        :  1. Active region is defined as [RemoveCell ... BoxSize-RemoveCell]
//                2. Useful only for non-periodic boundaries
//                   --> For removing particles lying too close to the simulation boundaries where
//                       the potential extrapolation may lead to large errors
//
// Parameter   :  x/y/z : Input coordinates
//
// Return      :  true/false : inside/outside the active region
//-------------------------------------------------------------------------------------------------------
bool Par_WithinActiveRegion( const real x, const real y, const real z )
{

// RemoveCell has the unit of root-level cell
   const double RemoveZone = amr->Par->RemoveCell*amr->dh[0];

// assuming OPT__BC_FLU[2*d] == OPT__BC_FLU[2*d+1] for d=0-2 when adoping periodic BC
   if (  OPT__BC_FLU[0] != BC_FLU_PERIODIC  &&  ( x < RemoveZone || x > amr->BoxSize[0]-RemoveZone )  )  return false;
   if (  OPT__BC_FLU[2] != BC_FLU_PERIODIC  &&  ( y < RemoveZone || y > amr->BoxSize[1]-RemoveZone )  )  return false;
   if (  OPT__BC_FLU[4] != BC_FLU_PERIODIC  &&  ( z < RemoveZone || z > amr->BoxSize[2]-RemoveZone )  )  return false;

   return true;

} // FUCNTION Par_WithinActiveRegion



#endif // #ifdef PARTICLE
