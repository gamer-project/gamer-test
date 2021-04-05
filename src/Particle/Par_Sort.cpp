#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SortByPos
// Description :  Sort particles by their position
//
// Note        :  1. Sorting by velocity may be necessary for STAR_FORMATION, where the new star particles
//                   created at different time but the same position may still have the same position for a
//                   while if velocity*dt is on the order of round-off errors
//                   --> Not supported yet since we may not have the velocity information (e.g., when adopting
//                       UseInputMassPos in Par_MassAssignment())
//                2. Currently IdxTable has the type "int" instead of "long" since Mis_Heapsort() doesn't
//                   support "long" yet...
//                3. Invoked by Par_MassAssignment() and FB_AdvanceDt()
//
// Parameter   :  NPar     : Number of particles
//                PosX/Y/Z : Particle position
//                IdxTable : Index table to be returned
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
void Par_SortByPos( const long NPar, const real *PosX, const real *PosY, const real *PosZ, int *IdxTable )
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

} // FUNCTION : Par_SortByPos



#endif // #ifdef PARTICLE
