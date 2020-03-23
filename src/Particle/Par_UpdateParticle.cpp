#include "GAMER.h"

#ifdef PARTICLE

//#ifndef GRAVITY
//#   error : ERROR : GRAVITY is not defined !!
//#endif

#ifdef GRAVITY
#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateParticle
// Description :  Update particle position and velocity at the target level
//
// Note        :  1. Does not work with "P5_Gradient"
//                2. Does not take into account the "periodic B.C." when updating particles
//                   --> After update, the particle position may lie outside the simulation box
//                   --> It will be corrected by Par_PassParticle2Sibling()
//                   --> However, periodic B.C. is taken into account when calculating gravity if OPT__BC_POT == BC_POT_PERIODIC
//                3. For the K-D-K scheme, this function performs either prediction (K-D) or correction (last K) operation
//                   --> Use the input parameter "UpdateStep" to control
//                4. For the Euler scheme, this function completes the full update
//                   --> UpdateStep==PAR_UPSTEP_CORR is meaningless
//                5. For KDK, particle time is set to -0.5*dt after the K-D operation to indicate that they require
//                   velocity correction (the last K operation) later. Otherwise particles just cross from fine to coarse
//                   grids cannot be distinguished from those already in the coarse grid, while we only want to apply
//                   velocity correction to the former. After the velocity correction, particle time is set to TimeNew
//                   For Euler, particle time is set to TimeNew in the PAR_UPSTEP_PRED step
//                6. Particle time may not be synchronized (so different particles may have different time).
//                   --> For example, particles just cross from coarse (lv) to fine (lv+1) grids may have time greater than
//                       other particles at lv+1. Also, particles just cross from fine (lv) to coarse (lv-1) grids may have
//                       time less than other particles at lv-1.
//                   --> Only update particles with time < TimeNew
//                7. If STORE_PAR_ACC is on, we can store the acceleration of each particle (option "StoreAcc")
//                   and also use the stored acceleration for advancing particles (option "UseStoredAcc"). Note that
//                   it does NOT work with the Euler integration.
//                8. For PAR_UPSTEP_ACC_ONLY, we only calculate and store the acceleration of each particle
//                   --> Particle position, velocity, and time are not modified at all
//                   --> Use "TimeNew" to determine the target time
//                   --> StoreAcc must be on, and UseStoredAcc must be off
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach (also used by PAR_UPSTEP_ACC_ONLY)
//                TimeOld      : Physical time before update
//                UpdateStep   : PAR_UPSTEP_PRED/PAR_UPSTEP_CORR/PAR_UPSTEP_ACC_ONLY
//                               (CORR is for PAR_INTEG_KDK only, and ACC_ONLY is for STORE_PAR_ACC only)
//                StoreAcc     : Store the acceleration of each particle (must work with STORE_PAR_ACC)
//                UseStoredAcc : Use the acceleration of each particle stored previously for advancing particles
//                               (must work with STORE_PAR_ACC)
//-------------------------------------------------------------------------------------------------------
void Par_UpdateParticle( const int lv, const double TimeNew, const double TimeOld, const ParUpStep_t UpdateStep,
                         const bool StoreAcc, const bool UseStoredAcc )
{

#  ifndef GRAVITY
   return;
#  else

   const ParInterp_t IntScheme    = amr->Par->Interp;
   const bool   UsePot            = ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT );
   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const double dh                = amr->dh[lv];
   const double _dh               = 1.0/dh;
   const double PrepPotTime       = ( UpdateStep==PAR_UPSTEP_CORR || UpdateStep==PAR_UPSTEP_ACC_ONLY ) ? TimeNew : TimeOld;

   const int  GraGhost_Par        = 1;                      // always set to 1 for particles (P5 gradient is not supported yet)
   const int  ParGhost            = amr->Par->GhostSize;
   const int  PotGhost            = GraGhost_Par + ParGhost;
   const int  PotSize             = PS1 + 2*PotGhost;
   const int  AccSize             = PS1 + 2*ParGhost;
   const real Const_8             = (real)8.0;
// const real GraConst            = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // but P5 is NOT supported yet
   const real GraConst            = ( false                ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // but P5 is NOT supported yet

   real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
#  ifdef STORE_PAR_ACC
   real *ParAcc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
#  endif
   real *ParTime   = amr->Par->Time;
   real *ParType   = amr->Par->Type;

// determine PotSg for STORE_POT_GHOST
#  ifdef STORE_POT_GHOST
   int  PotSg;
   real PotWeighting0, PotWeighting1;
   bool PotIntTime = false;

   if ( UsePot  &&  amr->Par->ImproveAcc )
   {
      if      (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][   amr->PotSg[lv] ], NULL, false )  )
         PotSg =   amr->PotSg[lv];

      else if (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][ 1-amr->PotSg[lv] ], NULL, false )  )
         PotSg = 1-amr->PotSg[lv];

      else
      {
         if ( OPT__DT_LEVEL == DT_LEVEL_SHARED )
         Aux_Error( ERROR_INFO, "Cannot determine PotSg for OPT__DT_LEVEL == DT_LEVEL_SHARED (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                    lv, PrepPotTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
         else if ( UpdateStep == PAR_UPSTEP_PRED )
         Aux_Error( ERROR_INFO, "Potential interpolation should be unnecessary (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                    lv, PrepPotTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );

         PotIntTime    = true;
         PotWeighting0 = ( +amr->PotSgTime[lv][1] - PrepPotTime ) / ( amr->PotSgTime[lv][1] - amr->PotSgTime[lv][0] );
         PotWeighting1 = ( -amr->PotSgTime[lv][0] + PrepPotTime ) / ( amr->PotSgTime[lv][1] - amr->PotSgTime[lv][0] );
         PotSg         = -1;   // useless
      }
   }
#  endif


// check
#  ifdef COMOVING
#  error : ERROR : does not support COMOVING because time-step has not been converted to comoving !!
#  endif

#  ifdef DEBUG_PARTICLE
// Par->ImproveAcc only works particle interpolation schemes with ParGhost == 1 (CIC & TSC)
   if ( amr->Par->ImproveAcc  &&  PotGhost != GRA_GHOST_SIZE )
      Aux_Error( ERROR_INFO, "PotGhost (%d) != GRA_GHOST_SIZE (%d) for amr->Par->ImproveAcc !!\n",
                 PotGhost, GRA_GHOST_SIZE );

// Par->ImproveAcc must work with STORE_POT_GHOST
#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "amr->Par->ImproveAcc must work with STORE_POT_GHOST !!\n" );
#  endif

   if ( amr->Par->Integ == PAR_INTEG_EULER  &&  UpdateStep == PAR_UPSTEP_CORR )
      Aux_Error( ERROR_INFO, "UpdateStep == PAR_UPSTEP_CORR is meaningless for PAR_INTEG_EULER !!\n" );

   if ( UpdateStep != PAR_UPSTEP_PRED  &&  UpdateStep != PAR_UPSTEP_CORR  &&  UpdateStep != PAR_UPSTEP_ACC_ONLY )
      Aux_Error( ERROR_INFO, "unsupported UpdateStep (%d) !!\n", UpdateStep );

#  ifdef STORE_PAR_ACC
   if ( StoreAcc  &&  amr->Par->Integ == PAR_INTEG_KDK  &&  UpdateStep == PAR_UPSTEP_PRED )
      Aux_Message( stderr, "WARNING : \"StoreAcc\" for PAR_INTEG_KDK and PAR_UPSTEP_PRED !?\n" );

   if ( StoreAcc  &&  amr->Par->Integ == PAR_INTEG_EULER  &&  UpdateStep == PAR_UPSTEP_CORR )
      Aux_Message( stderr, "WARNING : \"StoreAcc\" for PAR_INTEG_EULER and PAR_UPSTEP_CORR !?\n" );

   if ( UseStoredAcc  &&  amr->Par->Integ == PAR_INTEG_EULER )
      Aux_Message( stderr, "WARNING : \"UseStoredAcc\" for PAR_INTEG_EULER !?\n" );

   if ( UseStoredAcc  &&  amr->Par->Integ == PAR_INTEG_KDK  &&  UpdateStep == PAR_UPSTEP_CORR )
      Aux_Message( stderr, "WARNING : \"UseStoredAcc\" for PAR_INTEG_KDK and PAR_UPSTEP_CORR !?\n" );

   if ( UpdateStep == PAR_UPSTEP_ACC_ONLY )
   {
      if ( !StoreAcc )     Aux_Error( ERROR_INFO, "PAR_UPSTEP_ACC_ONLY must work with StoreAcc !!\n" );
      if ( UseStoredAcc )  Aux_Error( ERROR_INFO, "PAR_UPSTEP_ACC_ONLY does NOT work with UseStoredAcc !!\n" );
   }
#  else
   if ( UpdateStep == PAR_UPSTEP_ACC_ONLY )
      Aux_Error( ERROR_INFO, "\"PAR_UPSTEP_ACC_ONLY\" does not work when STORE_PAR_ACC if off !!\n" );

   if ( StoreAcc )      Aux_Error( ERROR_INFO, "\"StoreAcc\" does not work when STORE_PAR_ACC if off !!\n" );

   if ( UseStoredAcc )  Aux_Error( ERROR_INFO, "\"UseStoredAcc\" does not work when STORE_PAR_ACC if off !!\n" );
#  endif
#  endif // #ifdef DEBUG_PARTICLE


// OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   real *Pot = new real [ 8*CUBE(PotSize) ];    // 8: number of patches per patch group
   real *Acc = new real [ 3*CUBE(AccSize) ];    // 3: three dimension

   real (*Pot3D)[PotSize][PotSize][PotSize] = ( real (*)[PotSize][PotSize][PotSize] )Pot;
   real (*Acc3D)[AccSize][AccSize][AccSize] = ( real (*)[AccSize][AccSize][AccSize] )Acc;

   bool GotYou;
   long ParID;
   real Acc_Temp[3], dt, dt_half;


// loop over all **real** patch groups
#  pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with target particles
//    --> use patch group as the calculation unit since Prepare_PatchData() only work with patch group
//    --> disadvantage: some patches may not have particles ... (they will be skipped later)
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         int NActivePar = 0;
         for (int i=1; i<PAR_NTYPE; i++) 
            NActivePar += amr->patch[0][lv][PID]->NParType[i];

         if ( NActivePar > 0 )
         {
            if ( UpdateStep == PAR_UPSTEP_CORR )
            {
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParTime[ParID] < (real)0.0 && ParType[ParID] != PTYPE_TRACER )
                  {
                     GotYou = true;
                     break;
                  }
               }
            }

            else // UpdateStep == PAR_UPSTEP_PRED  ||  UpdateStep == PAR_UPSTEP_ACC_ONLY
               GotYou = true;
         }

         if ( GotYou )  break;
      } // for (int PID=PID0; PID<PID0+8; PID++)


//    nothing to do if there are no target particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the potential data for the patch group with particles (need NSIDE_26 for ParGhost>0 )
      if ( !UseStoredAcc  &&  UsePot )
      {
#        ifdef STORE_POT_GHOST
         if ( amr->Par->ImproveAcc )
         {
            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
            {
               int NActivePar = 0;
               for (int i=1; i<PAR_NTYPE; i++) 
                  NActivePar += amr->patch[0][lv][PID]->NParType[i];
               if ( NActivePar == 0 )  continue;   // skip patches with no massive particles

//             temporal interpolation is required for correcting the velocity of particles just crossing
//             from fine to coarse grids
               if ( PotIntTime )
               {
                  for (int k=0; k<PotSize; k++)
                  for (int j=0; j<PotSize; j++)
                  for (int i=0; i<PotSize; i++)
                     Pot3D[P][k][j][i] =   PotWeighting0*amr->patch[0][lv][PID]->pot_ext[k][j][i]
                                         + PotWeighting1*amr->patch[1][lv][PID]->pot_ext[k][j][i];
               }

               else
                  memcpy( Pot3D[P], amr->patch[PotSg][lv][PID]->pot_ext, CUBE(PotSize)*sizeof(real) );
            }
         } // if ( amr->Par->ImproveAcc )

         else
#        endif // #ifdef STORE_POT_GHOST
            Prepare_PatchData( lv, PrepPotTime, Pot, NULL, PotGhost, 1, &PID0, _POTE, _NONE,
                               OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                               OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );
      } // if ( !UseStoredAcc  &&  UsePot )


      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
         int NActivePar = 0;
         for (int i=1; i<PAR_NTYPE; i++) 
            NActivePar += amr->patch[0][lv][PID]->NParType[i];

         if ( NActivePar == 0 )  continue;   // skip patches with no massive particles

         if ( !UseStoredAcc )
         {
            double PhyCorner_ExtAcc[3], x, y, z;

            for (int d=0; d<3; d++)
               PhyCorner_ExtAcc[d] = amr->patch[0][lv][PID]->EdgeL[d] + (0.5-ParGhost)*dh;

//          3. calculate acceleration on cells
            for (int k=GraGhost_Par, kk=0; k<PotSize-GraGhost_Par; k++, kk++)  {  z = PhyCorner_ExtAcc[2] + (double)kk*dh;
            for (int j=GraGhost_Par, jj=0; j<PotSize-GraGhost_Par; j++, jj++)  {  y = PhyCorner_ExtAcc[1] + (double)jj*dh;
            for (int i=GraGhost_Par, ii=0; i<PotSize-GraGhost_Par; i++, ii++)  {  x = PhyCorner_ExtAcc[0] + (double)ii*dh;

               Acc_Temp[0] = (real)0.0;
               Acc_Temp[1] = (real)0.0;
               Acc_Temp[2] = (real)0.0;

//             3.1 external acceleration
               if ( OPT__EXT_ACC )
                  CPUExtAcc_Ptr( Acc_Temp, x, y, z, PrepPotTime, ExtAcc_AuxArray );

//             3.2 self-gravity and external potential
               if ( UsePot )
               {
//                OPT__GRA_P5_GRADIENT is not supported yet
//                if ( OPT__GRA_P5_GRADIENT )
                  if ( false )
                  {
                     Acc_Temp[0] += GraConst * ( -         Pot3D[P][k  ][j  ][i+2] +         Pot3D[P][k  ][j  ][i-2]
                                                 + Const_8*Pot3D[P][k  ][j  ][i+1] - Const_8*Pot3D[P][k  ][j  ][i-1] );
                     Acc_Temp[1] += GraConst * ( -         Pot3D[P][k  ][j+2][i  ] +         Pot3D[P][k  ][j-2][i  ]
                                                 + Const_8*Pot3D[P][k  ][j+1][i  ] - Const_8*Pot3D[P][k  ][j-1][i  ] );
                     Acc_Temp[2] += GraConst * ( -         Pot3D[P][k+2][j  ][i  ] +         Pot3D[P][k-2][j  ][i  ]
                                                 + Const_8*Pot3D[P][k+1][j  ][i  ] - Const_8*Pot3D[P][k-1][j  ][i  ] );
                  }

                  else
                  {
                     Acc_Temp[0] += GraConst * ( Pot3D[P][k  ][j  ][i+1] - Pot3D[P][k  ][j  ][i-1] );
                     Acc_Temp[1] += GraConst * ( Pot3D[P][k  ][j+1][i  ] - Pot3D[P][k  ][j-1][i  ] );
                     Acc_Temp[2] += GraConst * ( Pot3D[P][k+1][j  ][i  ] - Pot3D[P][k-1][j  ][i  ] );
                  }
               } // if ( UsePot )

               for (int d=0; d<3; d++)    Acc3D[d][kk][jj][ii] = Acc_Temp[d];

            }}} // i,j,k
         } // if ( !UseStoredAcc )


         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip tracer particles
            if ( ParType[ParID] == PTYPE_TRACER )
               continue;

//          determine time-step and skip particles with zero or negative time-step
            if ( UpdateStep == PAR_UPSTEP_PRED )
            {
//             it's crucial to first calculate dt here and skip particles with dt <= (real)0.0 (including the equal sign)
//             since later on we select particles with negative particle time (which has been set to -dt), with equal sign
//             excluded, for the velocity correction
               dt      = (real)TimeNew - ParTime[ParID];
               dt_half = (real)0.5*dt;

               if ( dt <= (real)0.0 )  continue;
            }

            else if ( UpdateStep == PAR_UPSTEP_CORR )
            {
//             during the prediction step, we store particle time as -0.5*dt (which must be < 0.0) to indicate that
//             these particles require velocity correction
               dt_half = -ParTime[ParID];

               if ( dt_half <= (real)0.0 )   continue;
            }

            else if ( UpdateStep == PAR_UPSTEP_ACC_ONLY )
            {
               dt      = NULL_REAL;    // useless
               dt_half = NULL_REAL;    // useless
            }


//          4. calculate acceleration at the particle position
            switch ( IntScheme ) {

//          4.1 NGP
            case ( PAR_INTERP_NGP ):
            {
               int idx[3];

//             calculate the nearest grid index
               if ( !UseStoredAcc )
               for (int d=0; d<3; d++)
               {
                  idx[d] = int( ( ParPos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh );

//                prevent from round-off errors (especially for NGP and TSC)
                  if ( idx[d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                                d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idx[d] );
#                    endif

                     idx[d] = 0;
                  }

                  else if ( idx[d] >= AccSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                                   d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idx[d] );
#                    endif

                     idx[d] = AccSize - 1;
                  }
               } // for (int d=0; d<3; d++)

//             calculate acceleration
               for (int d=0; d<3; d++)
               {
#                 ifdef STORE_PAR_ACC
                  if ( UseStoredAcc )
                     Acc_Temp[d] = ParAcc[d][ParID];
                  else
#                 endif
                     Acc_Temp[d] = Acc3D[d][ idx[2] ][ idx[1] ][ idx[0] ];

#                 ifdef STORE_PAR_ACC
                  if ( StoreAcc )   ParAcc[d][ParID] = Acc_Temp[d];
#                 endif
               }
            } // PAR_INTERP_NGP
            break;


//          4.2 CIC
            case ( PAR_INTERP_CIC ):
            {
               int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
               double dr      [3];     // distance to the center of the left cell
               double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

               if ( !UseStoredAcc )
               for (int d=0; d<3; d++)
               {
//                calculate the array index of the left and right cells
                  dr      [d] = ( ParPos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost - 0.5;
                  idxLR[0][d] = int( dr[d] );
                  idxLR[1][d] = idxLR[0][d] + 1;

//                prevent from round-off errors
//                (CIC should be clear off this issue unless round-off erros are comparable to dh)
                  if ( idxLR[0][d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                     idxLR[0][d] = 0;
                     idxLR[1][d] = 1;
                  }

                  else if ( idxLR[1][d] >= AccSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                     idxLR[0][d] = AccSize - 2;
                     idxLR[1][d] = AccSize - 1;
                  }

//                get the weighting of the nearby 8 cells
                  dr     [d] -= (double)idxLR[0][d];
                  Frac[0][d]  = 1.0 - dr[d];
                  Frac[1][d]  =       dr[d];
               } // for (int d=0; d<3; d++)

//             calculate acceleration
               for (int d=0; d<3; d++)
               {
#                 ifdef STORE_PAR_ACC
                  if ( UseStoredAcc )
                     Acc_Temp[d] = ParAcc[d][ParID];
                  else
#                 endif
                  {
                     Acc_Temp[d] = (real)0.0;

                     for (int k=0; k<2; k++)
                     for (int j=0; j<2; j++)
                     for (int i=0; i<2; i++)
                     Acc_Temp[d] += Acc3D[d][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  }

#                 ifdef STORE_PAR_ACC
                  if ( StoreAcc )   ParAcc[d][ParID] = Acc_Temp[d];
#                 endif
               }
            } // PAR_INTERP_CIC
            break;


//          4.3 TSC
            case ( PAR_INTERP_TSC ):
            {
               int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
               double dr       [3];    // distance to the left edge of the central cell
               double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

               if ( !UseStoredAcc )
               for (int d=0; d<3; d++)
               {
//                calculate the array index of the left, central, and right cells
                  dr       [d] = ( ParPos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost;
                  idxLCR[1][d] = int( dr[d] );
                  idxLCR[0][d] = idxLCR[1][d] - 1;
                  idxLCR[2][d] = idxLCR[1][d] + 1;

//                prevent from round-off errors (especially for NGP and TSC)
                  if ( idxLCR[0][d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                     idxLCR[0][d] = 0;
                     idxLCR[1][d] = 1;
                     idxLCR[2][d] = 2;
                  }

                  else if ( idxLCR[2][d] >= AccSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( ParPos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the acc array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                d, ParPos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                     idxLCR[0][d] = AccSize - 3;
                     idxLCR[1][d] = AccSize - 2;
                     idxLCR[2][d] = AccSize - 1;
                  }

//                get the weighting of the nearby 27 cells
                  dr     [d] -= (double)idxLCR[1][d];
                  Frac[0][d]  = 0.5*SQR( 1.0 - dr[d] );
                  Frac[1][d]  = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
                  Frac[2][d]  = 0.5*SQR( dr[d] );
               } // for (int d=0; d<3; d++)

//             calculate acceleration
               for (int d=0; d<3; d++)
               {
#                 ifdef STORE_PAR_ACC
                  if ( UseStoredAcc )
                     Acc_Temp[d] = ParAcc[d][ParID];
                  else
#                 endif
                  {
                     Acc_Temp[d] = (real)0.0;

                     for (int k=0; k<3; k++)
                     for (int j=0; j<3; j++)
                     for (int i=0; i<3; i++)
                     Acc_Temp[d] += Acc3D[d][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  }

#                 ifdef STORE_PAR_ACC
                  if ( StoreAcc )   ParAcc[d][ParID] = Acc_Temp[d];
#                 endif
               }
            } // PAR_INTERP_TSC
            break;


            default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
            } // switch ( IntScheme )


//          5. update particles
//          5.0 nothing to do if we only want to store particle acceleration
            if ( UpdateStep == PAR_UPSTEP_ACC_ONLY )     continue;


//          5.1 Euler method
            else if ( amr->Par->Integ == PAR_INTEG_EULER )
            {
               for (int d=0; d<3; d++)
               {
                  ParPos[d][ParID] += ParVel  [d][ParID]*dt;   // update position first
                  ParVel[d][ParID] += Acc_Temp[d]       *dt;
               }

               ParTime[ParID] = TimeNew;
            }


//          5.2 KDK scheme
            else if ( amr->Par->Integ == PAR_INTEG_KDK )
            {
//             5.2.1 KDK prediction
               if ( UpdateStep == PAR_UPSTEP_PRED )
               {
                  for (int d=0; d<3; d++)
                  {
                     ParVel[d][ParID] += Acc_Temp[d]       *dt_half; // predict velocity for 0.5*dt
                     ParPos[d][ParID] += ParVel  [d][ParID]*dt;      // update position by the half-step velocity for a full dt
                  }

                  ParTime[ParID] = -dt_half;   // negative --> indicating that it requires velocity correction
               }

//             5.2.2 KDK correction for velocity
               else // UpdateStep == PAR_UPSTEP_CORR
               {
                  for (int d=0; d<3; d++)
                     ParVel[d][ParID] += Acc_Temp[d]*dt_half;     // correct velocity for 0.5*dt

                  ParTime[ParID] = TimeNew;
               }
            } // amr->Par->Integ
         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)`
      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

// 6. free memory
   delete [] Pot;
   delete [] Acc;

   } // end of OpenMP parallel region

#endif // #ifndef GRAVITY

} // FUNCTION : Par_UpdateParticle



#endif // #ifdef PARTICLE
