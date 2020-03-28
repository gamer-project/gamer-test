#include "GAMER.h"

#ifdef PARTICLE

extern double Advect_Vel[3];
extern int    Advect_NPar[3];

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_AdvectTracers
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

void Par_Init_ByFunction_AdvectTracers( const long NPar_ThisRank, const long NPar_AllRank,
                                        real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                        real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                        real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   const real delta_p[3] = { amr->BoxSize[0]/(Advect_NPar[0]+1), 
                             amr->BoxSize[1]/(Advect_NPar[1]+1),
                             amr->BoxSize[2]/(Advect_NPar[2]+1) };

   const real delta_box = amr->BoxSize[0] / MPI_NRank;

   const long NParZ_Local = Advect_NPar[2] / MPI_NRank;
   const long NPar_ToAssign = Advect_NPar[0]*Advect_NPar[1]*NParZ_Local;

   if ( NPar_ToAssign != NPar_ThisRank )
      Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                 NPar_ToAssign, NPar_ThisRank );

   MPI_Barrier( MPI_COMM_WORLD );

// Generate particles within the uniform sphere of high energy

   for (long kk=0; kk<NParZ_Local;    kk++) 
   for (long jj=0; jj<Advect_NPar[1]; jj++) 
   for (long ii=0; ii<Advect_NPar[0]; ii++) 
   {

      const long p = IDX321( ii, jj, kk, Advect_NPar[0], Advect_NPar[1] );

      ParMass[p] = 0.0;

      ParPosX[p] = (ii+1)*delta_p[0];
      ParPosY[p] = (jj+1)*delta_p[1];
      ParPosZ[p] = (kk+1)*delta_p[2] + MPI_Rank*delta_box;

      ParVelX[p] = Advect_Vel[0];
      ParVelY[p] = Advect_Vel[1];
      ParVelZ[p] = Advect_Vel[2];

//    synchronize all particles to the physical time at the base level
      ParTime[p] = Time[0];

      ParType[p] = PTYPE_TRACER;

   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_AdvectTracers

#endif // #ifdef PARTICLE





