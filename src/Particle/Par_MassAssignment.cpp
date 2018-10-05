#include "GAMER.h"

#ifdef PARTICLE

static bool WithinRho( const int idx[], const int RhoSize );
static bool FarAwayParticle( real ParPosX, real ParPosY, real ParPosZ, const bool Periodic[], const real PeriodicSize_Phy[],
                             const real EdgeL[], const real EdgeR[] );
#ifdef BITWISE_REPRODUCIBILITY
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, int *IdxTable );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_MassAssignment
// Description :  Deposit particle mass onto grid
//
// Note        :  1. Three different schemes are supported
//                   (1) NGP : Nearest-Grid-Point
//                   (2) CIC : Could-In-Cell
//                   (3) TSC : Triangular-Shaped-Cloud
//                2. The deposited density field will be stored in Rho[]
//                   --> This array will be initialized as zero only if "InitZero=true"
//                3. Particles having no contribution to Rho[] (which has the range EdgeL[d] <= r[d] < EdgeL[d]+RhoSize*dh )
//                   will be ignored
//                4. Particle position will be predicted to the target physical time if PredictPos is on
//                   --> But they will NOT be stored back to the global Pos[] array
//                   --> Also remember to skip particles waiting for velocity correction since they have time
//                       temporarily set to -dt (moreover, they should already be synchronized with TargetTime)
//                5. This function does NOT work with periodic boundary condition when the root grid has only
//                   two patches (i.e., one patch group) along any direction (i.e., NX0_TOT[0/1/2] = 2*PATCH_SIZE)
//                   --> It is because for that extreme case
//                       --> we can have RhoSize > PeriodicSize[] due to the ghost zones
//                       --> need to consider both the target particles and their image particles when assigning
//                           mass to Rho[] (in other words, each target particle may contribute to more than one cell
//                           even with the NGP scheme), which is not considered here!
//                   --> This is the reason for the check "if ( Periodic[d]  &&  RhoSize > PeriodicSize[d] ) ..."
//                6. For bitwise reproducibility, particles are sorted by their position before mass deposition
//                   --> Also refer to the note of the routine SortParticle[]
//
// Parameter   :  ParList         : List of target particle IDs
//                NPar            : Number of particles
//                IntScheme       : Particle interpolation scheme
//                Rho             : Array to store the output density field (assumed to be a cubic array)
//                RhoSize         : Size of Rho[] along each direction
//                EdgeL           : Left edge of Rho[]
//                dh              : cell size of Rho[]
//                PredictPos      : true --> predict particle position to TargetTime
//                TargetTime      : Target time for predicting the particle position
//                InitZero        : True --> initialize Rho[] as zero
//                Periodic        : True --> apply periodic boundary condition to the target direction
//                PeriodicSize    : Number of cells in the periodic box (in the unit of dh)
//                UnitDens        : Assign unit density to each particle regardless of the real particle mass
//                                  and cell size
//                                  --> Useful for counting the number of particles on each cell (together
//                                      with IntScheme == PAR_INTERP_NGP), which can be used as a refinement
//                                      criterion (i.e., OPT__FLAG_NPAR_CELL)
//                CheckFarAway    : True --> check whether the input particles are far away from the given density array
//                                       --> If true, don't calculate their mass assignment cell indices at all
//                                       --> This may improve performance when some of the input particles have no
//                                           contribution at all to Rho[]. However, it may also introduce additional
//                                           overhead if most input particles do have contribution to Rho[].
//                UseInputMassPos : Use the input array InputMassPos[] to obtain particle mass and position
//                                  --> Used by LOAD_BALANCE, where particle position and mass may be stored in the
//                                      ParMassPos_Copy[] array of each patch
//                                  --> ParList[] becomes useless and must be set to NULL
//                                  --> Does not work with PredictPos since we don't have the information of particle
//                                      time and velocity
//                InputMassPos    : Particle mass and position arrays used by UseInputMassPos
//
// Return      :  Rho
//-------------------------------------------------------------------------------------------------------
void Par_MassAssignment( const long *ParList, const long NPar, const ParInterp_t IntScheme, real *Rho,
                         const int RhoSize, const double *EdgeL, const double dh, const bool PredictPos,
                         const double TargetTime, const bool InitZero, const bool Periodic[], const int PeriodicSize[3],
                         const bool UnitDens, const bool CheckFarAway, const bool UseInputMassPos, real **InputMassPos )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( Rho == NULL )   Aux_Error( ERROR_INFO, "Rho == NULL !!\n" );

   if ( EdgeL == NULL )    Aux_Error( ERROR_INFO, "EdgeL == NULL !!\n" );

   if ( PredictPos  &&  TargetTime < 0.0 )   Aux_Error( ERROR_INFO, "TargetTime = %14.7e < 0.0 !!\n", TargetTime );

   for (int d=0; d<3; d++)
      if ( Periodic[d]  &&  RhoSize > PeriodicSize[d] )
            Aux_Error( ERROR_INFO, "RhoSize (%d) > PeriodicSize[%d] (%d) !!\n", RhoSize, d, PeriodicSize[d] );

   if ( UseInputMassPos )
   {
      if ( ParList != NULL )  Aux_Error( ERROR_INFO, "ParList != NULL for UseInputMassPos !!\n" );
      if ( PredictPos )    Aux_Error( ERROR_INFO, "PredictPos does NOT work with UseInputMassPos !!\n" );
      if ( NPar > 0  &&  InputMassPos == NULL )
         Aux_Error( ERROR_INFO, "InputMassPos == NULL for UseInputMassPos (NPar = %ld) !!\n", NPar );

#     ifndef LOAD_BALANCE
      Aux_Message( stderr, "WARNING : are you sure you want to use UseInputMassPos when LOAD_BALANCE is off !?\n" );
#     endif
   }

   else
   {
      if ( NPar > 0  &&  ParList == NULL )   Aux_Error( ERROR_INFO, "ParList == NULL for NPar = %ld !!\n", NPar );
      if ( InputMassPos != NULL )   Aux_Error( ERROR_INFO, "InputMassPos != NULL when UseInputMassPos is off !!\n" );
   }
#  endif // #ifdef DEBUG_PARTICLE


// 1. initialization
   if ( InitZero )
      for (int t=0; t<CUBE(RhoSize); t++)    Rho[t] = (real)0.0;

   if ( NPar == 0 )  return;


// 2. set up attribute arrays, copy particle position since they might be modified during the position prediction
   real *Mass   = NULL;
   real *Pos[3] = { NULL, NULL, NULL };
   real *PType   = NULL;
   long  ParID, Idx;

   if ( UseInputMassPos )
   {
      Mass   = InputMassPos[PAR_MASS];
      Pos[0] = InputMassPos[PAR_POSX];
      Pos[1] = InputMassPos[PAR_POSY];
      Pos[2] = InputMassPos[PAR_POSZ];
      PType  = InputMassPos[PAR_TYPE];
   }

   else
   {
      Mass = new real [NPar];

      for (int d=0; d<3; d++)    Pos[d] = new real [NPar];

      PType = new real [NPar];

      for (long p=0; p<NPar; p++)
      {
         ParID = ParList[p];

         Mass  [p] = amr->Par->Mass[ParID];
         Pos[0][p] = amr->Par->PosX[ParID];
         Pos[1][p] = amr->Par->PosY[ParID];
         Pos[2][p] = amr->Par->PosZ[ParID];
         PType [p] = amr->Par->Type[ParID];
      }
   }


// 3. predict particle position
   if ( PredictPos )    Par_PredictPos( NPar, ParList, Pos[0], Pos[1], Pos[2], TargetTime );


// 3-1/2: sort particles by their position to fix the order of mass assignment
//        --> necessary for achieving bitwise reproducibility
#  ifdef BITWISE_REPRODUCIBILITY
   int *Sort_IdxTable = new int [NPar];   // it will fail if "long" is actually required for NPar

   SortParticle( NPar, Pos[0], Pos[1], Pos[2], Sort_IdxTable );
#  endif


// 4. deposit particle mass
   const double _dh       = 1.0 / dh;
   const double _dh3      = CUBE(_dh);
   const double Ghost_Phy = amr->Par->GhostSize*dh;

   real (*Rho3D)[RhoSize][RhoSize] = ( real (*)[RhoSize][RhoSize] )Rho;

   int  idx[3];      // array index for Rho
   real ParDens;     // mass density of the cloud
   real EdgeWithGhostL[3], EdgeWithGhostR[3], PeriodicSize_Phy[3];

   for (int d=0; d<3; d++)
   {
      EdgeWithGhostL  [d] = real( EdgeL[d] - Ghost_Phy );
      EdgeWithGhostR  [d] = real( EdgeL[d] + Ghost_Phy + RhoSize*dh );

      if ( Periodic[d] )
      PeriodicSize_Phy[d] = real( PeriodicSize[d]*dh );
   }


   switch ( IntScheme )
   {
//    4.1 NGP
      case ( PAR_INTERP_NGP ):
      {
         for (long p=0; p<NPar; p++)
         {
#           ifdef BITWISE_REPRODUCIBILITY
            Idx = Sort_IdxTable[p];
#           else
            Idx = p;
#           endif

//          4.1.0 ignore tracer particles
            if ( PType[Idx] == PTYPE_TRACER ) 
               continue;

//          4.1.1 discard particles far away from the target region
            if (  CheckFarAway  &&  FarAwayParticle( Pos[0][Idx], Pos[1][Idx], Pos[2][Idx],
                                                     Periodic, PeriodicSize_Phy, EdgeWithGhostL, EdgeWithGhostR )  )
               continue;

//          4.1.2 calculate the nearest grid index
            for (int d=0; d<3; d++)
            {
               idx[d] = (int)FLOOR( ( Pos[d][Idx] - EdgeL[d] )*_dh );

//             periodicity
               if ( Periodic[d] )
               {
                  idx[d] = ( idx[d] + PeriodicSize[d] ) % PeriodicSize[d];

#                 ifdef DEBUG_PARTICLE
                  if ( idx[d] < 0  ||  idx[d] >= PeriodicSize[d] )
                     Aux_Error( ERROR_INFO, "incorrect idx[%d] = %d (PeriodicSize = %d) !!\n",
                                d, idx[d], PeriodicSize[d] );
#                 endif
               }
            }

//          4.1.3 assign mass if within Rho[]
//          check inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( Mass[Idx] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", Idx, Mass[Idx] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = Mass[Idx]*_dh3;

            if (  WithinRho( idx, RhoSize )  )
               Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens;
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_NGP
      break;


//    4.2 CIC
      case ( PAR_INTERP_CIC ):
      {
         int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
         double dr      [3];     // distance to the center of the left cell
         double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

         for (long p=0; p<NPar; p++)
         {
#           ifdef BITWISE_REPRODUCIBILITY
            Idx = Sort_IdxTable[p];
#           else
            Idx = p;
#           endif

//          4.2.0 ignore tracer particles
            if ( PType[Idx] == PTYPE_TRACER ) 
               continue;

//          4.2.1 discard particles far away from the target region
            if (  CheckFarAway  &&  FarAwayParticle( Pos[0][Idx], Pos[1][Idx], Pos[2][Idx],
                                                     Periodic, PeriodicSize_Phy, EdgeWithGhostL, EdgeWithGhostR )  )
               continue;

            for (int d=0; d<3; d++)
            {
//             4.2.2 calculate the array index of the left and right cells
               dr      [d]  = ( Pos[d][Idx] - EdgeL[d] )*_dh - 0.5;
               idxLR[0][d]  = (int)FLOOR( dr[d] );
               idxLR[1][d]  = idxLR[0][d] + 1;
               dr      [d] -= (double)idxLR[0][d];

//             periodicity
               if ( Periodic[d] )
               {
                  for (int t=0; t<2; t++)
                  {
                     idxLR[t][d] = ( idxLR[t][d] + PeriodicSize[d] ) % PeriodicSize[d];

#                    ifdef DEBUG_PARTICLE
                     if ( idxLR[t][d] < 0  ||  idxLR[t][d] >= PeriodicSize[d] )
                        Aux_Error( ERROR_INFO, "incorrect idxLR[%d][%d] = %d (PeriodicSize = %d) !!\n",
                                   t, d, idxLR[t][d], PeriodicSize[d] );
#                    endif
                  }
               }

//             4.2.3 get the weighting of the nearby 8 cells
               Frac[0][d] = 1.0 - dr[d];
               Frac[1][d] =       dr[d];
            } // for (int d=0; d<3; d++)

//          4.2.4 assign mass if within Rho[]
//          check inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( Mass[Idx] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", Idx, Mass[Idx] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = Mass[Idx]*_dh3;

            for (int k=0; k<2; k++) {  idx[2] = idxLR[k][2];
            for (int j=0; j<2; j++) {  idx[1] = idxLR[j][1];
            for (int i=0; i<2; i++) {  idx[0] = idxLR[i][0];

               if (  WithinRho( idx, RhoSize )  )
                  Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens*Frac[i][0]*Frac[j][1]*Frac[k][2];

            }}}
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_CIC
      break;


//    4.3 TSC
      case ( PAR_INTERP_TSC ):
      {
         int    idxLCR[3][3];    // array index of the left (idxLCR[0][d]), central (idxLCR[1][d]) and right (idxLCR[2][d]) cells
         double dr       [3];    // distance to the left edge of the central cell
         double Frac  [3][3];    // weighting of the left (Frac[0][d]), central (Frac[1][d]) and right (Frac[2][d]) cells

         for (long p=0; p<NPar; p++)
         {
#           ifdef BITWISE_REPRODUCIBILITY
            Idx = Sort_IdxTable[p];
#           else
            Idx = p;
#           endif

//          4.3.0 ignore tracer particles
            if ( PType[Idx] == PTYPE_TRACER ) 
               continue;

//          4.3.1 discard particles far away from the target region
            if (  CheckFarAway  &&  FarAwayParticle( Pos[0][Idx], Pos[1][Idx], Pos[2][Idx],
                                                     Periodic, PeriodicSize_Phy, EdgeWithGhostL, EdgeWithGhostR )  )
               continue;

            for (int d=0; d<3; d++)
            {
//             4.3.2 calculate the array index of the left, central, and right cells
               dr       [d]  = ( Pos[d][Idx] - EdgeL[d] )*_dh;
               idxLCR[1][d]  = (int)FLOOR( dr[d] );
               idxLCR[0][d]  = idxLCR[1][d] - 1;
               idxLCR[2][d]  = idxLCR[1][d] + 1;
               dr       [d] -= (double)idxLCR[1][d];

//             periodicity
               if ( Periodic[d] )
               {
                  for (int t=0; t<3; t++)
                  {
                     idxLCR[t][d] = ( idxLCR[t][d] + PeriodicSize[d] ) % PeriodicSize[d];

#                    ifdef DEBUG_PARTICLE
                     if ( idxLCR[t][d] < 0  ||  idxLCR[t][d] >= PeriodicSize[d] )
                        Aux_Error( ERROR_INFO, "incorrect idxLCR[%d][%d] = %d (PeriodicSize = %d) !!\n",
                                   t, d, idxLCR[t][d], PeriodicSize[d] );
#                    endif
                  }
               }

//             4.3.3 get the weighting of the nearby 27 cells
               Frac[0][d] = 0.5*SQR( 1.0 - dr[d] );
               Frac[1][d] = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
               Frac[2][d] = 0.5*SQR( dr[d] );
            } // for (int d=0; d<3; d++)

//          4.3.4 assign mass if within Rho[]
//          check inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( Mass[Idx] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", Idx, Mass[Idx] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = Mass[Idx]*_dh3;

            for (int k=0; k<3; k++) {  idx[2] = idxLCR[k][2];
            for (int j=0; j<3; j++) {  idx[1] = idxLCR[j][1];
            for (int i=0; i<3; i++) {  idx[0] = idxLCR[i][0];

               if (  WithinRho( idx, RhoSize )  )
                  Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens*Frac[i][0]*Frac[j][1]*Frac[k][2];
            }}}
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_TSC
      break;

      default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
   } // switch ( IntScheme )


// 5. free memory
   if ( !UseInputMassPos )
   {
      delete [] Mass;
      for (int d=0; d<3; d++)    delete [] Pos[d];
   }

#  ifdef BITWISE_REPRODUCIBILITY
   delete [] Sort_IdxTable;
#  endif

} // FUNCTION : Par_MassAssignment



//-------------------------------------------------------------------------------------------------------
// Function    :  WithinRho
// Description :  Check whether the input cell lies within the Rho array
//
// Note        :  This function does NOT consider periodicity
//
// Parameter   :  idx     : Cell indices
//                RhoSize : Size of the Rho array
//
// Return      :  true --> within the target region
//-------------------------------------------------------------------------------------------------------
bool WithinRho( const int idx[], const int RhoSize )
{

   if ( idx[0] < 0  ||  idx[0] >= RhoSize  ||
        idx[1] < 0  ||  idx[1] >= RhoSize  ||
        idx[2] < 0  ||  idx[2] >= RhoSize    )
      return false;

   else
      return true;

} // FUNCTION : WithinRho



//-------------------------------------------------------------------------------------------------------
// Function    :  FarAwayParticle
// Description :  Check whether the input particle is far away from the target density array and thus
//                cannot have any contribution
//
// Note        :  1. Periodic BC is taken care of by considering the position of image particles
//                2. Particles pass this check are guaranteed to have contribution to the density array.
//                   However, some cells with mass deposited from these particles may still lie outside
//                   the density array. Therefore, we must still call WithinRho() to check furthre.
//
// Parameter   :  ParPosX/Y/Z      : Particle position
//                Periodic         : True --> apply periodic boundary condition to the target direction
//                PeriodicSize_Phy : Size of the periodic box
//                EdgeL/R          : Left and right edge of the density array (including ghost zones outside
//                                   the density array)
//
// Return      :  (true / false) <--> particle (has no / does have) contribution to the give density array
//-------------------------------------------------------------------------------------------------------
bool FarAwayParticle( real ParPosX, real ParPosY, real ParPosZ, const bool Periodic[], const real PeriodicSize_Phy[],
                      const real EdgeL[], const real EdgeR[] )
{

// x
   if ( Periodic[0] )
   {
      if      ( ParPosX <  EdgeL[0] )  ParPosX += PeriodicSize_Phy[0];
      else if ( ParPosX >= EdgeR[0] )  ParPosX -= PeriodicSize_Phy[0];
   }

   if ( ParPosX < EdgeL[0]  ||  ParPosX >= EdgeR[0] )    return true;


// y
   if ( Periodic[1] )
   {
      if      ( ParPosY <  EdgeL[1] )  ParPosY += PeriodicSize_Phy[1];
      else if ( ParPosY >= EdgeR[1] )  ParPosY -= PeriodicSize_Phy[1];
   }

   if ( ParPosY < EdgeL[1]  ||  ParPosY >= EdgeR[1] )    return true;


// z
   if ( Periodic[2] )
   {
      if      ( ParPosZ <  EdgeL[2] )  ParPosZ += PeriodicSize_Phy[2];
      else if ( ParPosZ >= EdgeR[2] )  ParPosZ -= PeriodicSize_Phy[2];
   }

   if ( ParPosZ < EdgeL[2]  ||  ParPosZ >= EdgeR[2] )    return true;


   return false;

} // FUNCTION : FarAwayParticle



#ifdef BITWISE_REPRODUCIBILITY
//-------------------------------------------------------------------------------------------------------
// Function    :  SortParticle
// Description :  Sort particles by their position
//
// Note        :  1. Sorting by velocity may be necessary for STAR_FORMATION, where the new star particles
//                   created at different time but the same position may still have the same position for a
//                   while if velocity*dt is on the order of round-off errors
//                   --> Not supported yet since we may not have the velocity information (e.g., when InputMassPos is adopted)
//                2. Currently IdxTable has the type "int" instead of "long" since Mis_Heapsort() doesn't
//                   support "long" yet...
//
// Parameter   :  NPar     : Number of particles
//                PosX/Y/Z : Particle position
//                IdxTable : Index table to be returned
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
void SortParticle( const long NPar, const real *PosX, const real *PosY, const real *PosZ, int *IdxTable )
{

   long NParSameX, NParSameY;

// 0. back up the PosX array since we don't want to modify it during the sorting
   real *PosX_Sorted = new real [NPar];

   memcpy( PosX_Sorted, PosX, sizeof(real)*NPar );


// 1. sort by x
   Mis_Heapsort( NPar, PosX_Sorted, IdxTable );


// 2. sort by y
   for (long x=0; x<NPar-1; x++)
   {
      NParSameX = 1;    // 1 --> itself

//    find all particles with the same x coordinate
      while ( x+NParSameX<NPar  &&  PosX_Sorted[x] == PosX_Sorted[x+NParSameX] )    NParSameX ++;

      if ( NParSameX > 1 )
      {
         int  *SortByX_IdxTable = new int  [NParSameX];  // it will fail if "long" is actually required
         int  *SortByY_IdxTable = new int  [NParSameX];  // it will fail if "long" is actually required
         real *PosY_Sorted      = new real [NParSameX];

         for (long y=0; y<NParSameX; y++)
         {
            SortByX_IdxTable[y] = IdxTable[ x + y ];
            PosY_Sorted     [y] = PosY[ SortByX_IdxTable[y] ];
         }

         Mis_Heapsort( NParSameX, PosY_Sorted, SortByY_IdxTable );

         for (long y=0; y<NParSameX; y++)    IdxTable[ x + y ] = SortByX_IdxTable[ SortByY_IdxTable[y] ];


//       3. sort by z
         for (long y=0; y<NParSameX-1; y++)
         {
            NParSameY = 1;    // 1 --> itself

//          find all particles with the same x and y coordinates
            while ( y+NParSameY<NParSameX  &&  PosY_Sorted[y] == PosY_Sorted[y+NParSameY] )  NParSameY ++;

            if ( NParSameY > 1 )
            {
               int  *SortByZ_IdxTable = new int  [NParSameY];  // it will fail if "long" is actually required
               real *PosZ_Sorted      = new real [NParSameY];

               for (long z=0; z<NParSameY; z++)
               {
                  SortByY_IdxTable[z] = IdxTable[ x + y + z ];
                  PosZ_Sorted     [z] = PosZ[ SortByY_IdxTable[z] ];
               }

               Mis_Heapsort( NParSameY, PosZ_Sorted, SortByZ_IdxTable );

               for (long z=0; z<NParSameY; z++)    IdxTable[ x + y + z ] = SortByY_IdxTable[ SortByZ_IdxTable[z] ];

               delete [] SortByZ_IdxTable;
               delete [] PosZ_Sorted;

               y += NParSameY - 1;
            } // if ( NParSameY > 1 )
         } // for (long y=0; y<NParSameX-1; y++)

         delete [] SortByX_IdxTable;
         delete [] SortByY_IdxTable;
         delete [] PosY_Sorted;

         x += NParSameX - 1;
      } // if ( NParSameX > 1 )
   } // for (long x=0; x<NPar-1; x++)

   delete [] PosX_Sorted;

} // FUNCTION : SortParticle
#endif // #ifdef BITWISE_REPRODUCIBILITY



#endif // #ifdef PARTICLE
