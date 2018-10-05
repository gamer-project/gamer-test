#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_GetConservedQuantity
// Description :  Calculate conserved quantities for particles
//                --> Mass, momentum, kinematic energy and potential energy
//
// Note        :  1. Use "call by reference" to return results
//                2. Return Ep=0.0 if GRAVITY is off
//                3. External potential is included, but external acceleration is NOT included
//                4. Particles may NOT be full synchronized when calculating their energies
//                   --> But it is found to have a minimal effect since most particles are synchronized
//                   --> We don't want to synchronize particles in this routine (by calling Par_Synchronize)
//                       since particles may move outside the current patch which requires additional workload
//                5. Results obtained from serial (with SERIAL on) and parallel codes (with LOAD_BALANCE) on
//                   can be slightly different due to round-off errors
//                   --> We do not correct it even in the debug mode because it should have no impact on the
//                       actual simulation variables
//                   --> Different number of OpenMP threads can also results in different results again
//                       due to round-off errors
//
// Parameter   :  Mass_Total     : Total particle mass to be returned
//                MomX/Y/Z_Total : Total momentum to be returned
//                Ek_Total       : Total kinematic energy to be returned
//                Ep_Total       : Total potential energy to be returned
//
// Return      :  Mass_Total, MomX/Y/Z_Total, Ek_Total, Ep_Total
//-------------------------------------------------------------------------------------------------------
void Par_Aux_GetConservedQuantity( double &Mass_Total, double &MomX_Total, double &MomY_Total, double &MomZ_Total,
                                   double &Ek_Total, double &Ep_Total )
{

// 1. mass, momentum, and kinematic energy
   double Mass_ThisRank=0.0, MomX_ThisRank=0.0, MomY_ThisRank=0.0, MomZ_ThisRank=0.0, Ek_ThisRank=0.0;
   double Send[5], Recv[5];

// use static schedule to give the same reduction results everytime
#  pragma omp parallel for schedule( static ) reduction( +:Mass_ThisRank, MomX_ThisRank, MomY_ThisRank, MomZ_ThisRank, Ek_ThisRank )
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
//    skip inactive, massless, and tracer particles
      if ( amr->Par->Mass[p] > 0.0 || amr->Par->Type[p] == PTYPE_TRACER)
      {
         Mass_ThisRank += amr->Par->Mass[p];
         MomX_ThisRank += amr->Par->Mass[p]*amr->Par->VelX[p];
         MomY_ThisRank += amr->Par->Mass[p]*amr->Par->VelY[p];
         MomZ_ThisRank += amr->Par->Mass[p]*amr->Par->VelZ[p];
         Ek_ThisRank   += 0.5*amr->Par->Mass[p]*( SQR(amr->Par->VelX[p]) + SQR(amr->Par->VelY[p]) + SQR(amr->Par->VelZ[p]) );
      }
   }

   Send[0] = Mass_ThisRank;
   Send[1] = MomX_ThisRank;
   Send[2] = MomY_ThisRank;
   Send[3] = MomZ_ThisRank;
   Send[4] = Ek_ThisRank;

   MPI_Reduce( Send, Recv, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

   Mass_Total = Recv[0];
   MomX_Total = Recv[1];
   MomY_Total = Recv[2];
   MomZ_Total = Recv[3];
   Ek_Total   = Recv[4];


// 2. potential energy
#  ifdef GRAVITY
   const ParInterp_t IntScheme  = amr->Par->Interp;
   const bool UsePot            = ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT );
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const real MinTemp_No        = -1.0;

   const int  PotGhost          = amr->Par->GhostSize;
   const int  PotSize           = PS1 + 2*PotGhost;

   const real *Pos[3]           = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *Mass             = amr->Par->Mass;
   const real *PType            = amr->Par->Type;
 
   double Ep_ThisRank = 0.0;
   double PrepPotTime, dh, _dh, Ep_Coeff;
   int    PotSg;

// set potential energy to zero when enabling both OPT__SELF_GRAVITY and OPT__EXT_POT
// since the potential energy obtained here would be wrong anyway
// --> to avoid possible misinterpretation
   if      (  OPT__SELF_GRAVITY  &&  !OPT__EXT_POT )  Ep_Coeff = 0.5;
   else if ( !OPT__SELF_GRAVITY  &&   OPT__EXT_POT )  Ep_Coeff = 1.0;
   else                                               Ep_Coeff = 0.0;


// check
// Par->ImproveAcc only supports particle interpolation schemes with ParGhost == 1 (CIC & TSC)
   if ( amr->Par->ImproveAcc  &&  PotGhost != GRA_GHOST_SIZE-1 )
      Aux_Error( ERROR_INFO, "PotGhost (%d) != GRA_GHOST_SIZE-1 (%d) for amr->Par->ImproveAcc !!\n",
                 PotGhost, GRA_GHOST_SIZE-1 );

// Par->ImproveAcc must work with STORE_POT_GHOST
#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "amr->Par->ImproveAcc must work with STORE_POT_GHOST !!\n" );
#  endif


// nothing to do if not using potential
   if ( !UsePot )    Ep_Total = 0.0;

// loop over particles in all leaf patches to compute potential energy
   else
   for (int lv=0; lv<NLEVEL; lv++)
   {
      dh          = amr->dh[lv];
      _dh         = 1.0/dh;
      PrepPotTime = Time[lv];

//    determine PotSg for STORE_POT_GHOST
#     ifdef STORE_POT_GHOST
      if ( amr->Par->ImproveAcc )
      {
         if      (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][   amr->PotSg[lv] ], NULL, false )  )
            PotSg =   amr->PotSg[lv];

         else if (  Mis_CompareRealValue( PrepPotTime, amr->PotSgTime[lv][ 1-amr->PotSg[lv] ], NULL, false )  )
            PotSg = 1-amr->PotSg[lv];

         else
            Aux_Error( ERROR_INFO, "Cannot determine PotSg (lv %d, PrepTime %20.14e, SgTime[0] %20.14e, SgTime[1] %20.14e !!\n",
                       lv, PrepPotTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1] );
      }
#     endif


//    OpenMP parallel region
#     pragma omp parallel
      {

//    per-thread variables
      bool   GotYou;
      long   ParID;

      real *Pot = new real [ 8*CUBE(PotSize) ];    // 8: number of patches per patch group
      real (*Pot3D)[PotSize][PotSize][PotSize] = ( real (*)[PotSize][PotSize][PotSize] )Pot;


//    use static schedule to give the same reduction result everytime
#     pragma omp for schedule( static ) reduction( +:Ep_ThisRank )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
//       2-1. find the patch groups with particles
//       --> use patch group as the calculation unit since Prepare_PatchData() only work with patch group
//       --> some patches may not have particles ...
         GotYou = false;

         for (int PID=PID0; PID<PID0+8; PID++)
         {
            if ( amr->patch[0][lv][PID]->NPar > 0 )
            {
               GotYou = true;
               break;
            }
         }


//       nothing to do if there are no particles in the target patch group
         if ( !GotYou )    continue;


//       2-2. prepare the potential data for the patch group with particles (need NSIDE_26 for ParGhost>0 )
#        ifdef STORE_POT_GHOST
         if ( amr->Par->ImproveAcc )
         {
            const int didx = 1;  // assuming GRA_GHOST_SIZE - Pot_Size = 1

            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
            {
               if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

               for (int k=didx; k<GRA_NXT-didx; k++)
               for (int j=didx; j<GRA_NXT-didx; j++)
               for (int i=didx; i<GRA_NXT-didx; i++)
                  Pot3D[P][k-didx][j-didx][i-didx] = amr->patch[PotSg][lv][PID]->pot_ext[k][j][i];
            }
         }

         else
#        endif
            Prepare_PatchData( lv, PrepPotTime, Pot, NULL, PotGhost, 1, &PID0, _POTE, _NONE,
                               OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, (PotGhost==0)?NSIDE_00:NSIDE_26, IntPhase_No,
                               OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );


//       2-3. calculate potential energy
         for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
         {
            if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;   // skip patches with no particles

            switch ( IntScheme )
            {
//             2-3-1 NGP
               case ( PAR_INTERP_NGP ):
               {
                  int idx[3];

                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
                  {
                     ParID = amr->patch[0][lv][PID]->ParList[p];

                     if ( PType[ParID] == PTYPE_TRACER || Mass[ParID] < 0.0 )
                        continue;

//                   calculate the nearest grid index
                     for (int d=0; d<3; d++)
                     {
                        idx[d] = int( ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh );

//                      prevent from round-off errors (especially for NGP and TSC)
                        if ( idx[d] < 0 )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idx[d] );
#                          endif

                           idx[d] = 0;
                        }

                        else if ( idx[d] >= PotSize )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                              Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                                         d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idx[d] );
#                          endif

                           idx[d] = PotSize - 1;
                        }
                     } // for (int d=0; d<3; d++)

//                   get potential energy
                     Ep_ThisRank += Ep_Coeff*Mass[ParID]*Pot3D[P][ idx[2] ][ idx[1] ][ idx[0] ];
                  } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               } // PAR_INTERP_NGP
               break;


//             2-3-2. CIC
               case ( PAR_INTERP_CIC ):
               {
                  int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
                  double dr      [3];     // distance to the center of the left cell
                  double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
                  {
                     ParID = amr->patch[0][lv][PID]->ParList[p];
                     
                     if ( PType[ParID] == PTYPE_TRACER || Mass[ParID] < 0.0 )
                        continue;

                     for (int d=0; d<3; d++)
                     {
//                      calculate the array index of the left and right cells
                        dr      [d]  = ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + PotGhost - 0.5;
                        idxLR[0][d]  = int( dr[d] );
                        idxLR[1][d]  = idxLR[0][d] + 1;
                        dr      [d] -= (double)idxLR[0][d];

//                      prevent from round-off errors
//                      (CIC should be clear off this issue unless round-off erros are comparable to dh)
                        if ( idxLR[0][d] < 0 )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLR[0][d], idxLR[1][d] );
#                          endif

                           idxLR[0][d] = 0;
                           idxLR[1][d] = 1;
                        }

                        else if ( idxLR[1][d] >= PotSize )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLR[0][d], idxLR[1][d] );
#                          endif

                           idxLR[0][d] = PotSize - 2;
                           idxLR[1][d] = PotSize - 1;
                        }

//                      get the weighting of the nearby 8 cells
                        Frac[0][d] = 1.0 - dr[d];
                        Frac[1][d] =       dr[d];
                     } // for (int d=0; d<3; d++)

//                   get potential energy
                     for (int k=0; k<2; k++)
                     for (int j=0; j<2; j++)
                     for (int i=0; i<2; i++)
                     Ep_ThisRank += Ep_Coeff*Mass[ParID]*Pot3D[P][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
                                            *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               } // PAR_INTERP_CIC
               break;


//             2-3-3. TSC
               case ( PAR_INTERP_TSC ):
               {
                  int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
                  double dr       [3];    // distance to the left edge of the central cell
                  double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
                  {
                     ParID = amr->patch[0][lv][PID]->ParList[p];

                     if ( PType[ParID] == PTYPE_TRACER || Mass[ParID] < 0.0 )
                        continue;

                     for (int d=0; d<3; d++)
                     {
//                      calculate the array index of the left, central, and right cells
                        dr       [d]  = ( Pos[d][ParID] - amr->patch[0][lv][PID]->EdgeL[d] )*_dh + PotGhost;
                        idxLCR[1][d]  = int( dr[d] );
                        idxLCR[0][d]  = idxLCR[1][d] - 1;
                        idxLCR[2][d]  = idxLCR[1][d] + 1;
                        dr       [d] -= (double)idxLCR[1][d];

//                      prevent from round-off errors (especially for NGP and TSC)
                        if ( idxLCR[0][d] < 0 )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeL[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#                          endif

                           idxLCR[0][d] = 0;
                           idxLCR[1][d] = 1;
                           idxLCR[2][d] = 2;
                        }

                        else if ( idxLCR[2][d] >= PotSize )
                        {
#                          ifdef DEBUG_PARTICLE
                           if (  ! Mis_CompareRealValue( Pos[d][ParID], (real)amr->patch[0][lv][PID]->EdgeR[d], NULL, false )  )
                           Aux_Error( ERROR_INFO, "index outside the pot array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                                      d, Pos[d][ParID], amr->patch[0][lv][PID]->EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#                          endif

                           idxLCR[0][d] = PotSize - 3;
                           idxLCR[1][d] = PotSize - 2;
                           idxLCR[2][d] = PotSize - 1;
                        }

//                      get the weighting of the nearby 27 cells
                        Frac[0][d] = 0.5*SQR( 1.0 - dr[d] );
                        Frac[1][d] = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
                        Frac[2][d] = 0.5*SQR( dr[d] );
                     } // for (int d=0; d<3; d++)

//                   get potential energy
                     for (int k=0; k<3; k++)
                     for (int j=0; j<3; j++)
                     for (int i=0; i<3; i++)
                     Ep_ThisRank += Ep_Coeff*Mass[ParID]*Pot3D[P][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
                                            *Frac[i][0]*Frac[j][1]*Frac[k][2];
                  } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
               } // PAR_INTERP_TSC
               break;

               default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
            } // switch ( IntScheme )
         } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

//    2-4. free memory
      delete [] Pot;

      } // end of OpenMP parallel region
   } // for (int lv=0; lv<NLEVEL; lv++)


// 2.5. gather from all ranks
   MPI_Reduce( &Ep_ThisRank, &Ep_Total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

#  else // #ifdef GRAVITY

   Ep_Total = 0.0;
#  endif // #ifdef GRAVITY ... else ...

} // FUNCTION : Par_Aux_GetConservedQuantity



#endif // #ifdef PARTICLE
