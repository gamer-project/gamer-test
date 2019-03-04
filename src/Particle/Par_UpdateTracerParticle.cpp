#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateTracerParticle
// Description :  Update tracer particle position and velocity at the target level
//
// Note        :  1. Does not take into account the "periodic B.C." when updating particles
//                   --> After update, the particle position may lie outside the simulation box
//                   --> It will be corrected by Par_PassParticle2Sibling()
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
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach (also used by PAR_UPSTEP_ACC_ONLY)
//                TimeOld      : Physical time before update
//-------------------------------------------------------------------------------------------------------
void Par_UpdateTracerParticle( const int lv, const double TimeNew, const double TimeOld)
{

   const ParInterp_t IntScheme    = amr->Par->InterpTracer;
   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const double dh                = amr->dh[lv];
   const double _dh               = 1.0/dh;

   const int  VelGhost_Par        = 1;                      // always set to 1 for particles 
   const int  ParGhost            = amr->Par->GhostSize;
   const int  PotGhost            = VelGhost_Par + ParGhost;
   const int  VelSize             = PS1 + 2*ParGhost;

   real InterpParPos[3];

   real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   real *ParTime   = amr->Par->Time;
   real *ParType   = amr->Par->Type;



// check
#  ifdef COMOVING
#  error : ERROR : does not support COMOVING because time-step has not been converted to comoving !!
#  endif

// OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   real *VelX = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group
   real *VelY = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group
   real *VelZ = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group

   real (*VelX3D)[VelSize][VelSize][VelSize] = ( real (*)[VelSize][VelSize][VelSize] )VelX;
   real (*VelY3D)[VelSize][VelSize][VelSize] = ( real (*)[VelSize][VelSize][VelSize] )VelY;
   real (*VelZ3D)[VelSize][VelSize][VelSize] = ( real (*)[VelSize][VelSize][VelSize] )VelZ;

   bool   GotYou;
   long   ParID;
   real   Vel_Temp[3], dt, dt_half;
   double x, y, z;


// loop over all **real** patch groups
#  pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with target tracer particles
//    --> use patch group as the calculation unit since Prepare_PatchData() only work with patch group
//    --> disadvantage: some patches may not have tracer particles ... (they will be skipped later)
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         if ( amr->patch[0][lv][PID]->NParType[1] > 0 )
         {
            if ( UpdateStep == PAR_UPSTEP_CORR )
            {
               for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               {
                  ParID = amr->patch[0][lv][PID]->ParList[p];

                  if ( ParTime[ParID] < (real)0.0 && ParType[ParID] == PTYPE_TRACER )
                  {
                     GotYou = true;
                     break;
                  }
               }
            }

            else // UpdateStep == PAR_UPSTEP_PRED 
               GotYou = true;
         }

         if ( GotYou )  break;
      } // for (int PID=PID0; PID<PID0+8; PID++)


//    nothing to do if there are no target tracer particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the velocity data for the patch group with particles (need NSIDE_26 for ParGhost>0 )

      Prepare_PatchData( lv, TimeNew, Pot, PotGhost, 1, &PID0, _POTE,
                         OPT__GRA_INT_SCHEME, UNIT_PATCH, NSIDE_26, IntPhase_No, 
                         OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
         if ( amr->patch[0][lv][PID]->NParType[1] == 0 )  continue;   // skip patches with no tracer particles

         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip massive particles
            if ( ParType[ParID] == PTYPE_MASSIVE )
               continue;

//          determine time-step and skip particles with zero or negative time-step
            dt      = (real)TimeNew - ParTime[ParID];
            if ( dt <= (real)0.0 )  continue;
            dt_half = (real)0.5*dt;

            for (int d=0; d<3; d++) InterpParPos[d] = ParPos[d][ParID]+dt*ParVel[d][ParID];
            
//          4. calculate gas velocity at the particle position
            switch ( IntScheme ) {

//          4.1 NGP
            case ( PAR_INTERP_NGP ):
            {
               int idx[3];

//             calculate the nearest grid index
               for (int d=0; d<3; d++)
               {
                  idx[d] = int( ( InterpParPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh );

//                prevent from round-off errors (especially for NGP and TSC)
                  if ( idx[d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                                d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idx[d] );
#                    endif

                     idx[d] = 0;
                  }

                  else if ( idx[d] >= VelSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                        Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                                   d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idx[d] );
#                    endif

                     idx[d] = VelSize - 1;
                  }
               } // for (int d=0; d<3; d++)

//             calculate new particle velocity
               Vel_Temp[0] = VelX3D[P][ idx[2] ][ idx[1] ][ idx[0] ];
               Vel_Temp[1] = VelY3D[P][ idx[2] ][ idx[1] ][ idx[0] ];
               Vel_Temp[2] = VelZ3D[P][ idx[2] ][ idx[1] ][ idx[0] ];

            } // PAR_INTERP_NGP
            break;


//          4.2 CIC
            case ( PAR_INTERP_CIC ):
            {
               int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
               double dr      [3];     // distance to the center of the left cell
               double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

               for (int d=0; d<3; d++)
               {
//                calculate the array index of the left and right cells
                  dr      [d] = ( InterpParPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost - 0.5;
                  idxLR[0][d] = int( dr[d] );
                  idxLR[1][d] = idxLR[0][d] + 1;

//                prevent from round-off errors
//                (CIC should be clear off this issue unless round-off erros are comparable to dh)
                  if ( idxLR[0][d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                     idxLR[0][d] = 0;
                     idxLR[1][d] = 1;
                  }

                  else if ( idxLR[1][d] >= VelSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                     idxLR[0][d] = VelSize - 2;
                     idxLR[1][d] = VelSize - 1;
                  }

//                get the weighting of the nearby 8 cells
                  dr     [d] -= (double)idxLR[0][d];
                  Frac[0][d]  = 1.0 - dr[d];
                  Frac[1][d]  =       dr[d];
               } // for (int d=0; d<3; d++)

//             calculate velocity
               for (int d=0; d<3; d++) Vel_Temp[d] = (real)0.0;

               for (int k=0; k<2; k++) {
                  for (int j=0; j<2; j++) {
                     for (int i=0; i<2; i++) {
                        Vel_Temp[0] += VelX3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                        Vel_Temp[1] += VelY3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                        Vel_Temp[2] += VelZ3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                     }
                  }
               }

            } // PAR_INTERP_CIC
            break;


//          4.3 TSC
            case ( PAR_INTERP_TSC ):
            {
               int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
               double dr       [3];    // distance to the left edge of the central cell
               double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

               for (int d=0; d<3; d++)
               {
//                calculate the array index of the left, central, and right cells
                  dr       [d] = ( InterpParPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost;
                  idxLCR[1][d] = int( dr[d] );
                  idxLCR[0][d] = idxLCR[1][d] - 1;
                  idxLCR[2][d] = idxLCR[1][d] + 1;

//                prevent from round-off errors (especially for NGP and TSC)
                  if ( idxLCR[0][d] < 0 )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                     idxLCR[0][d] = 0;
                     idxLCR[1][d] = 1;
                     idxLCR[2][d] = 2;
                  }

                  else if ( idxLCR[2][d] >= VelSize )
                  {
#                    ifdef DEBUG_PARTICLE
                     if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                     Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                     idxLCR[0][d] = VelSize - 3;
                     idxLCR[1][d] = VelSize - 2;
                     idxLCR[2][d] = VelSize - 1;
                  }

//                get the weighting of the nearby 27 cells
                  dr     [d] -= (double)idxLCR[1][d];
                  Frac[0][d]  = 0.5*SQR( 1.0 - dr[d] );
                  Frac[1][d]  = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
                  Frac[2][d]  = 0.5*SQR( dr[d] );
               } // for (int d=0; d<3; d++)

//             calculate velocity
               for (int d=0; d<3; d++) Vel_Temp[d] = (real)0.0;

               for (int k=0; k<3; k++) {
                  for (int j=0; j<3; j++) {
                     for (int i=0; i<3; i++) {
                        Vel_Temp[0] += VelX3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                        Vel_Temp[1] += VelY3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                        Vel_Temp[2] += VelZ3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                   *Frac[i][0]*Frac[j][1]*Frac[k][2];
                     }
                  }
               }

            } // PAR_INTERP_TSC
            break;


            default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
            } // switch ( IntScheme )


//          5. update particles


//          5.0 Euler method
            if ( amr->Par->IntegTracer == PAR_INTEG_EULER )
            {
               for (int d=0; d<3; d++)
               {
                  ParPos[d][ParID] = InterpParPos[d];
                  ParVel[d][ParID] = Vel_Temp[d];
               }

            }

//          5.1 RK2 scheme
            else if ( amr->Par->IntegTracer == PAR_INTEG_RK2 )
            {

               for (int d=0; d<3; d++)
               {
                  ParPos[d][ParID] += dt_half*(ParVel[d][ParID]+Vel_Temp[d]);
                  ParVel[d][ParID] = Vel_Temp[d]; // This is wrong now, need to fix it
               }

            } // amr->Par->IntegTracer

            ParTime[ParID] = TimeNew;

         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)`
      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

// 6. free memory
   delete [] VelX;
   delete [] VelY;
   delete [] VelZ;

   } // end of OpenMP parallel region

} // FUNCTION : Par_UpdateTracerParticle

void interpolate_velocity(real InterpPos[3], real InterpVel[3], ) 
{

// 4. calculate gas velocity at the particle position
   switch ( IntScheme ) {

//    4.1 NGP
      case ( PAR_INTERP_NGP ):
      {
         int idx[3];

//       calculate the nearest grid index
         for (int d=0; d<3; d++)
         {
            idx[d] = int( ( InterpPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh );

//          prevent from round-off errors (especially for NGP and TSC)
            if ( idx[d] < 0 )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idx[d] );
#                    endif

                  idx[d] = 0;
            }

            else if ( idx[d] >= VelSize )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idx[d] );
#                    endif

                  idx[d] = VelSize - 1;
            }
            } // for (int d=0; d<3; d++)

//             calculate new particle velocity
            Vel_Temp[0] = VelX3D[P][ idx[2] ][ idx[1] ][ idx[0] ];
            Vel_Temp[1] = VelY3D[P][ idx[2] ][ idx[1] ][ idx[0] ];
            Vel_Temp[2] = VelZ3D[P][ idx[2] ][ idx[1] ][ idx[0] ];

      } // PAR_INTERP_NGP
      break;


//          4.2 CIC
      case ( PAR_INTERP_CIC ):
      {
            int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
            double dr      [3];     // distance to the center of the left cell
            double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

            for (int d=0; d<3; d++)
            {
//                calculate the array index of the left and right cells
            dr      [d] = ( InterpParPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost - 0.5;
            idxLR[0][d] = int( dr[d] );
            idxLR[1][d] = idxLR[0][d] + 1;

//                prevent from round-off errors
//                (CIC should be clear off this issue unless round-off erros are comparable to dh)
            if ( idxLR[0][d] < 0 )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                  idxLR[0][d] = 0;
                  idxLR[1][d] = 1;
            }

            else if ( idxLR[1][d] >= VelSize )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idxLR[0][d], idxLR[1][d] );
#                    endif

                  idxLR[0][d] = VelSize - 2;
                  idxLR[1][d] = VelSize - 1;
            }

//                get the weighting of the nearby 8 cells
            dr     [d] -= (double)idxLR[0][d];
            Frac[0][d]  = 1.0 - dr[d];
            Frac[1][d]  =       dr[d];
            } // for (int d=0; d<3; d++)

//             calculate velocity
            for (int d=0; d<3; d++) Vel_Temp[d] = (real)0.0;

            for (int k=0; k<2; k++) {
            for (int j=0; j<2; j++) {
                  for (int i=0; i<2; i++) {
                  Vel_Temp[0] += VelX3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  Vel_Temp[1] += VelY3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  Vel_Temp[2] += VelZ3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  }
            }
            }

      } // PAR_INTERP_CIC
      break;


//          4.3 TSC
      case ( PAR_INTERP_TSC ):
      {
            int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
            double dr       [3];    // distance to the left edge of the central cell
            double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

            for (int d=0; d<3; d++)
            {
//                calculate the array index of the left, central, and right cells
            dr       [d] = ( InterpParPos[d] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + ParGhost;
            idxLCR[1][d] = int( dr[d] );
            idxLCR[0][d] = idxLCR[1][d] - 1;
            idxLCR[2][d] = idxLCR[1][d] + 1;

//                prevent from round-off errors (especially for NGP and TSC)
            if ( idxLCR[0][d] < 0 )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                  idxLCR[0][d] = 0;
                  idxLCR[1][d] = 1;
                  idxLCR[2][d] = 2;
            }

            else if ( idxLCR[2][d] >= VelSize )
            {
#                    ifdef DEBUG_PARTICLE
                  if (  ! Mis_CompareRealValue( InterpParPos[d], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the vel array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                              d, InterpParPos[d], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                    endif

                  idxLCR[0][d] = VelSize - 3;
                  idxLCR[1][d] = VelSize - 2;
                  idxLCR[2][d] = VelSize - 1;
            }

//                get the weighting of the nearby 27 cells
            dr     [d] -= (double)idxLCR[1][d];
            Frac[0][d]  = 0.5*SQR( 1.0 - dr[d] );
            Frac[1][d]  = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
            Frac[2][d]  = 0.5*SQR( dr[d] );
            } // for (int d=0; d<3; d++)

//             calculate velocity
            for (int d=0; d<3; d++) Vel_Temp[d] = (real)0.0;

            for (int k=0; k<3; k++) {
            for (int j=0; j<3; j++) {
                  for (int i=0; i<3; i++) {
                  Vel_Temp[0] += VelX3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  Vel_Temp[1] += VelY3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  Vel_Temp[2] += VelZ3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                              *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  }
            }
            }

      } // PAR_INTERP_TSC
      break;


      default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      } // switch ( IntScheme )

} // FUNCTION : interpolate_velocity


#endif // #ifdef PARTICLE
