#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Father
// Description :  Pass particles from sons to father
//
// Note        :  1. After calling this function, son patches will have no particles (NPar == 0)
//                2. This function should always be called before deleting son patches within the same patch group
//                   --> invoked by "Refine & LB_Refine_AllocateNewPatch"
//
// Parameter   :  FaLv  : Father's refinement level
//                FaPID : Father's patch ID
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Father( const int FaLv, const int FaPID )
{

#  ifdef DEBUG_PARTICLE
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect FaLv = %d !!\n", FaLv );

   if ( FaPID < 0  ||  FaPID >= amr->num[FaLv] )
      Aux_Error( ERROR_INFO, "incorrect FaPID = %d (FaLv %d, NPatch %d) !!\n", FaPID, FaLv, amr->num[FaLv] );
#  endif


   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;


#  ifdef DEBUG_PARTICLE
   if ( SonPID0 < -1 )
      Aux_Error( ERROR_INFO, "This function does NOT work with sons living abroad (FaLv %d, FaPID %d, SonPID0 %d) !!\n",
                 FaLv, FaPID, SonPID0 );
#  endif


// nothing to do if father has no son
   if ( SonPID0 == -1 )    return;


// 1. get the total number of particles in all sons
   int NParSon = 0;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)  NParSon += amr->patch[0][SonLv][SonPID]->NPar;


// nothing to do if sons have no particles
   if ( NParSon == 0 )  return;


// 2. gather the particle lists from all sons
   long *ParListSon = new long [NParSon];
   int  t = 0;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
   for (int p=0; p<amr->patch[0][SonLv][SonPID]->NPar; p++)    ParListSon[ t++ ] = amr->patch[0][SonLv][SonPID]->ParList[p];


// 3. add particles to father
//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
   const real *PType = amr->Par->Type;
#  ifdef DEBUG_PARTICLE
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   amr->patch[0][FaLv][FaPID]->AddParticle( NParSon, ParListSon, &amr->Par->NPar_Lv[FaLv], PType, ParPos, amr->Par->NPar_AcPlusInac, __FUNCTION__ );
#  else
   amr->patch[0][FaLv][FaPID]->AddParticle( NParSon, ParListSon, &amr->Par->NPar_Lv[FaLv], PType );
#  endif


// 4. remove particles in all sons
//###NOTE : No OpenMP since RemoveParticle will modify amr->Par->NPar_Lv[]
   const bool RemoveAllParticle = true;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      amr->patch[0][SonLv][SonPID]->RemoveParticle( NULL_INT, NULL, 
         &amr->Par->NPar_Lv[SonLv], RemoveAllParticle, PType );


// free memory
   delete [] ParListSon;

} // FUNCTION : Par_PassParticle2Father



#endif // #ifdef PARTICLE
