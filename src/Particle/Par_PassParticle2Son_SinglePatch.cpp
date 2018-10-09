#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Son_SinglePatch
// Description :  Pass particles from father to sons
//
// Note        :  1. After calling this function, father patch will have no particles (NPar == 0)
//                   --> Particles should always reside in "leaf" patches
//                2. This function should always be called after new son patches are allocated
//                   --> Invoked by Init_Refine(), Refine(), LB_Refine_AllocateNewPatch(), and LB_Init_Refine()
//                3. This function can also be used even when son patches already have particles
//                   --> The case where particles just cross a coarse-fine boundary (from coarse to fine)
//                   --> Invoked by Par_PassParticle2Sibling_AllPatch(), which is invoked
//                       by EvolveLevel() after the velocity correction in KDK
//                4. Target father patch can be either real or buffer
//                5. Does NOT work with sons living abroad
//                   --> Use Par_PassParticle2Son_MultiPatch() for that
//
// Parameter   :  FaLv  : Father refinement level
//                FaPID : Father patch ID
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Son_SinglePatch( const int FaLv, const int FaPID )
{

   const int NPar    = amr->patch[0][FaLv][FaPID]->NPar;
   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;


#  ifdef DEBUG_PARTICLE
   if ( SonPID0 < -1 )
      Aux_Error( ERROR_INFO, "This function does NOT work with sons living abroad (FaLv %d, FaPID %d, SonPID0 %d) !!\n",
                 FaLv, FaPID, SonPID0 );
#  endif


// nothing to do if father has no son or no particles at home
   if ( SonPID0 < 0  ||  NPar == 0 )   return;


// 1. allocate the new particle list for each son
   long *NewListForSon[8];
   int   NNewForSon[8];
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      NewListForSon[LocalID] = new long [NPar];    // this is the maximum array size required
      NNewForSon   [LocalID] = 0;
   }


// 2. find the home patches at SonLv for all particles
   const int     Octant[2][2][2] = {  { {0,1},{2,4} }, { {3,6},{5,7} }  };    // LocalID of sons at different octants
   const double *FaCen           = amr->patch[0][SonLv][SonPID0+7]->EdgeL;    // central coordinates of FaPID
   const real   *ParPos[3]       = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real   *PType           = amr->Par->Type;
   long ParID;
   int  ijk[3], LocalID;

   for (int p=0; p<NPar; p++)
   {
      ParID = amr->patch[0][FaLv][FaPID]->ParList[p];

      for (int d=0; d<3; d++)    ijk[d] = ( ParPos[d][ParID] < FaCen[d] ) ? 0 : 1;

      LocalID = Octant[ ijk[2] ][ ijk[1] ][ ijk[0] ];
      NewListForSon[LocalID][ NNewForSon[LocalID] ++ ] = ParID;
   } // for (int p=0; p<NPar; p++)


// 3. pass particles to each son
   int SonPID;
   for (int LocalID=0; LocalID<8; LocalID++)
   {
      SonPID = SonPID0 + LocalID;

//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
#     ifdef DEBUG_PARTICLE
      amr->patch[0][SonLv][SonPID]->AddParticle( NNewForSon[LocalID], NewListForSon[LocalID], &amr->Par->NPar_Lv[SonLv],
                                                 PType, ParPos, amr->Par->NPar_AcPlusInac, __FUNCTION__ );
#     else
      amr->patch[0][SonLv][SonPID]->AddParticle( NNewForSon[LocalID], NewListForSon[LocalID], &amr->Par->NPar_Lv[SonLv], PType );
#     endif
   }


// 4. remove particles in the father patch
//###NOTE : No OpenMP since RemoveParticle will modify amr->Par->NPar_Lv[]
   const bool RemoveAllParticle = true;
   amr->patch[0][FaLv][FaPID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[FaLv], RemoveAllParticle );


// free memory
   for (int LocalID=0; LocalID<8; LocalID++)    delete [] NewListForSon[LocalID];

} // FUNCTION : Par_PassParticle2Son_SinglePatch



#endif // #ifdef PARTICLE
