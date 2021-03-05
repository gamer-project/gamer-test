#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_GetTimeStep_VelAcc
// Description :  Estimate the evolution time-step from the maximum particle velocity and acceleration
//                --> dt_vel = DT__PARVEL*dh/v_max, where v_max = max(vx,vy,vz,all_particles)
//                    dt_acc = DT__PARACC*(dh/a_max)^0.5, where a_max = max(accx,accy,accz,all_particles)
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Particle acceleration criterion is used only when DT__PARACC > 0.0
//
// Parameter   :  dt_vel : Evolution time-step estimated from the particle velocity
//                         --> Call-by-reference
//                dt_acc : Evolution time-step estimated from the particle acceleration
//                         --> Call-by-reference
//                lv     : Target refinement level
//
// Return      :  dt_vel, dt_acc
//-------------------------------------------------------------------------------------------------------
void Par_GetTimeStep_VelAcc( double &dt_vel, double &dt_acc, const int lv )
{

#  ifdef COMOVING
#  error : ERROR : COMOVING is not supported yet !!
#  endif


   const real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
#  if ( defined STORE_PAR_ACC && defined GRAVITY ) 
   const real *Acc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
#  else
   const real *Acc[3] = { NULL, NULL, NULL };
#  endif
#  ifdef MASSIVE_PARTICLES
   const bool  UseAcc = ( DT__PARACC > 0.0 );
#  else
   const bool  UseAcc = false;
#  endif
#  ifdef OPENMP
   const int   NT     = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int   NT     = 1;
#  endif

#  ifndef STORE_PAR_ACC
   if ( UseAcc )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) > 0.0 when STORE_PAR_ACC is off !!\n", DT__PARACC );
#  endif

   real  MaxVel, MaxAcc;
   int   TID;
   long  ParID;

   real *MaxVel_OMP = new real [NT];
   real *MaxAcc_OMP = new real [NT];


// get the maximum particle velocity and acceleration at the target level
   MaxVel = (real)0.0;  // don't assign negative values since we assume it to be positive-definite
   MaxAcc = (real)0.0;

   for (int t=0; t<NT; t++)
   {
      MaxVel_OMP[t] = (real)0.0;
      MaxAcc_OMP[t] = (real)0.0;
   }

//###NOTE: OpenMP may not improve performance here
#  pragma omp parallel private( TID, ParID )
   {
#     ifdef OPENMP
      TID = omp_get_thread_num();
#     else
      TID = 0;
#     endif

#     pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

         for (int d=0; d<3; d++)
         MaxVel_OMP[TID] = MAX( MaxVel_OMP[TID], FABS(Vel[d][ParID]) );

         if ( UseAcc )
         for (int d=0; d<3; d++)
         MaxAcc_OMP[TID] = MAX( MaxAcc_OMP[TID], FABS(Acc[d][ParID]) );
      }
   } // end of OpenMP parallel region

// compare the maximum velocity and acceleration evaluated by different OMP threads
   for (int t=0; t<NT; t++)
   {
      MaxVel = MAX( MaxVel_OMP[t], MaxVel );

      if ( UseAcc )
      MaxAcc = MAX( MaxAcc_OMP[t], MaxAcc );
   }

   delete [] MaxVel_OMP;
   delete [] MaxAcc_OMP;


// get the time-step in this rank
   double dt_vel_local, dt_acc_local;
   dt_vel_local =       amr->dh[lv] / MaxVel;

   if ( UseAcc )
   dt_acc_local = sqrt( amr->dh[lv] / MaxAcc );


// get the minimum time-step in all ranks
   MPI_Allreduce( &dt_vel_local, &dt_vel, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

   if ( UseAcc )
   MPI_Allreduce( &dt_acc_local, &dt_acc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( !Aux_IsFinite(dt_vel)  &&  amr->Par->NPar_Lv[lv] > 0  &&  MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : time-step estimation by particle velocity is incorrect (dt_vel = %13.7e) !!\n", dt_vel );
      Aux_Message( stderr, "          --> Likely all particles have zero velocity\n" );

      if ( DT__PARVEL_MAX < 0.0 )
      Aux_Message( stderr, "          --> You might want to set DT__PARVEL_MAX properly\n" );
   }

   if ( UseAcc  &&  !Aux_IsFinite(dt_acc)  &&  amr->Par->NPar_Lv[lv] )
      Aux_Error( ERROR_INFO, "time-step estimation by particle acceleration is incorrect (dt_acc = %13.7e) !!\n", dt_acc );


// multiply by the safety factor
   dt_vel *= DT__PARVEL;
   if ( DT__PARVEL_MAX >= 0.0 )  dt_vel = MIN( dt_vel, DT__PARVEL_MAX );

   if ( UseAcc )
   dt_acc *= DT__PARACC;

} // FUNCTION : Par_GetTimeStep_VelAcc



#endif // #ifdef PARTICLE
