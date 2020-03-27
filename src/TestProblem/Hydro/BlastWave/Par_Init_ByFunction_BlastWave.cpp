#include "GAMER.h"

#ifdef PARTICLE

extern RandomNumber_t *RNG;

extern double Blast_Center[3];
extern double Blast_Radius;


double BlastWave_RandomNumber(RandomNumber_t *RNG, const double Min, const double Max);


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Merger
// Description :  Initialize all particle attributes for the merging cluster test
//                --> Modified from "Par_Init_ByFile.cpp"
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
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
//                4. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
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

void Par_Init_ByFunction_BlastWave( const long NPar_ThisRank, const long NPar_AllRank,
                                    real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                    real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                    real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// Generate particles within the uniform sphere of high energy

   for (long p=0; p<NPar_ThisRank; p++)
   {

      ParMass[p] = 0.0;

      double phi = BlastWave_RandomNumber(RNG, 0.0, 2.0*M_PI);
      double theta = acos(BlastWave_RandomNumber(RNG, -1.0, 1.0));
      double r = BlastWave_RandomNumber(RNG, 0.0, Blast_Radius);

      ParPosX[p] = real(Blast_Center[0] + r*cos(phi)*sin(theta));
      ParPosY[p] = real(Blast_Center[1] + r*sin(phi)*sin(theta));
      ParPosZ[p] = real(Blast_Center[2] + r*cos(theta));

      ParVelX[p] = 0.0;
      ParVelY[p] = 0.0;
      ParVelZ[p] = 0.0;

//    synchronize all particles to the physical time at the base level
      ParTime[p] = Time[0];

      ParType[p] = PTYPE_TRACER;

   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_BlastWave



#endif // #ifdef PARTICLE





