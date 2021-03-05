#include "GAMER.h"

#ifdef PARTICLE

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Par_Init_ByFunction_Template( const long NPar_ThisRank, const long NPar_AllRank,
                                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );

// this function pointer must be set by a test problem initializer
void (*Par_Init_ByFunction_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                 real *ParType, real *AllAttribute[PAR_NATT_TOTAL] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Template
// Description :  Template of user-specified particle initializer
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr",
//                   which must be set by a test problem initializer
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Template( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)   ParTime[p] = Time[0];


// initialize the particle creation time by an arbitrary negative value since it is
// only used for star particles created during evolution and is useless during initialization
#  ifdef STAR_FORMATION
   for (int p=0; p<NPar_ThisRank; p++)    AllAttribute[Idx_ParCreTime][p] = -1.0;
#  endif


// set other particle attributes
// ============================================================================================================
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };

// example : randomly initialize
   /*
   const uint RSeed     = 2;                                         // random seed
   const real MassMin   = 1.0e-2;                                    // minimum value of particle mass
   const real MassMax   = 1.0;                                       // maximum value of particle mass
   const real PosMin[3] = { 0.0, 0.0, 0.0 };                         // minimum value of particle position
   const real PosMax[3] = { real( amr->BoxSize[0]*(1.0-1.0e-5) ),    // maximum value of particle position
                            real( amr->BoxSize[1]*(1.0-1.0e-5) ),
                            real( amr->BoxSize[2]*(1.0-1.0e-5) ) };
   const real VelMin[3] = { -1.0, -1.0, -1.0 };                      // minimum value of particle velocity
   const real VelMax[3] = { +1.0, +1.0, +1.0 };                      // maximum value of particle velocity

   srand( RSeed );

   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParMass[p] = ( (real)rand()/RAND_MAX )*( MassMax - MassMin ) + MassMin;

      for (int d=0; d<3; d++)
      {
         ParPos[d][p] = ( (real)rand()/RAND_MAX )*( PosMax[d] - PosMin[d] ) + PosMin[d];
         ParVel[d][p] = ( (real)rand()/RAND_MAX )*( VelMax[d] - VelMin[d] ) + VelMin[d];
      }
   }
   */
// ============================================================================================================


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Template



#endif // #ifdef PARTICLE
