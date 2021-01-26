#include "GAMER.h"

extern double (*Mis_GetTimeStep_User_Ptr)( const int lv, const double dTime_dt );




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep
// Description :  Estimate the evolution time-step (dt) and the physical time interval (dTime) at the target
//                refinement level
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the physical time interval (dTime) instead of the evolution time-step (dt)
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dTime, the physical time interval == "delta(scale_factor)"
//                       in the comoving coordinates, back to dt in EvolveLevel()
//                2. For OPT__DT_USER, the function pointer "Mis_GetTimeStep_User_Ptr" must be set by a
//                   test problem initializer
//
// Parameter   :  lv                : Target refinement level
//                dTime_SyncFaLv    : dt to synchronize lv and lv-1
//                                    --> Only used for OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE
//                AutoReduceDtCoeff : dt coefficient used by AUTO_REDUCE_DT
//                                    --> The final dt will be multiplied by this factor
//
// Return      :  dTime_min
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep( const int lv, const double dTime_SyncFaLv, const double AutoReduceDtCoeff )
{

   const char  FileName[] = "Record__TimeStep";
   static bool FirstTime  = true;
   const int   NdTimeMax  = 20;

   char  (*dTime_Name)[MAX_STRING] = new char   [NdTimeMax][MAX_STRING];
   double *dTime                   = new double [NdTimeMax];

   int NdTime = 0;


// -1. return immediately if the target level has no patches
// =============================================================================================================
   if ( NPatchTotal[lv] == 0 )   return HUGE_NUMBER;



// 0. estimate the relation between the evolution time-step (dt) and the physical time interval (dTime)
// =============================================================================================================
   double dTime_dt;  // first derivative of dTime over dt
#  ifdef COMOVING
   dTime_dt = pow(  OMEGA_M0*pow( Time[lv], 3.0 ) + (1.0-OMEGA_M0)*pow( Time[lv], 6.0 ),  0.5  );
#  else
   dTime_dt = 1.0;
#  endif



// 1.1 CRITERION ONE : fluid solver
// =============================================================================================================
#  if   ( MODEL == HYDRO )
   dTime[NdTime] = dTime_dt * dt_InvokeSolver( DT_FLU_SOLVER, lv );
   sprintf( dTime_Name[NdTime++], "%s", "Hydro_CFL" );

#  ifdef SRHD
   if ( DT_SPEED_OF_LIGHT ) dTime[NdTime] = ( (Step==0)?DT__FLUID_INIT:DT__FLUID ) * amr->dh[lv];
   else                     dTime[NdTime] = dt_InvokeSolver( DT_FLU_SOLVER, lv );
   sprintf( dTime_Name[NdTime++], "%s", "Hydro_CFL" );
#  endif

#  elif ( MODEL == ELBDM )
   dTime[NdTime] = dTime_dt * ELBDM_GetTimeStep_Fluid( lv );
   sprintf( dTime_Name[NdTime++], "%s", "ELBDM_CFL" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// 1.2 CRITERION TWO : gravitational acceleration
// =============================================================================================================
#  ifdef GRAVITY
#  if   ( MODEL == HYDRO )
   dTime[NdTime] = dTime_dt * dt_InvokeSolver( DT_GRA_SOLVER, lv );
   sprintf( dTime_Name[NdTime++], "%s", "Hydro_Acc" );

#  elif ( MODEL == ELBDM )
   dTime[NdTime] = dTime_dt * ELBDM_GetTimeStep_Gravity( lv  );
   sprintf( dTime_Name[NdTime++], "%s", "ELBDM_Pot" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
#  endif // #ifdef GRAVITY


// 1.3 CRITERION THREE : maximum allowed variation of the expansion factor
// =============================================================================================================
#  ifdef COMOVING
   dTime[NdTime] = DT__MAX_DELTA_A * Time[lv];
   sprintf( dTime_Name[NdTime++], "%s", "Delta_A" );
#  endif


// 1.4 CRITERION FOUR : match the time of the next data dump
// =============================================================================================================
// DumpByTime : true --> dump data according to the physical time
#  ifdef PARTICLE
   const bool DumpData   = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS ||
                             OPT__OUTPUT_PAR_TEXT );
#  else
   const bool DumpData   = ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS );
#  endif
   const bool DumpByTime = (  DumpData  &&  ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT || OPT__OUTPUT_MODE == OUTPUT_USE_TABLE )  )
                           ? true : false;

   if ( DumpByTime )
   {
      dTime[NdTime] = DumpTime - Time[lv];

      if ( dTime[NdTime] <= 0.0 )
      {
         Aux_Message( stderr, "********************************************************************************\n" );
         Aux_Message( stderr, "ERROR : dTime (%20.14e) <= 0.0, something is wrong !!\n", dTime[NdTime] );
         Aux_Message( stderr, "        (DumpTime %20.14e, Time %20.14e, lv %d)\n", DumpTime, Time[lv], lv );
         Aux_Message( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n",
                      MPI_Rank, __FILE__, __LINE__, __FUNCTION__ );
         Aux_Message( stderr, "********************************************************************************\n" );
         MPI_Exit();
      }

      sprintf( dTime_Name[NdTime++], "%s", "Data_Dump" );
   }


// 1.5 CRITERION FIVE : match the program end time
// =============================================================================================================
   dTime[NdTime] = END_T - Time[lv];

   if ( dTime[NdTime] <= 0.0 )
   {
      Aux_Message( stderr, "********************************************************************************\n" );
      Aux_Message( stderr, "ERROR : dTime (%20.14e) <= 0.0, something is wrong !!\n", dTime[NdTime] );
      Aux_Message( stderr, "        (END_T %20.14e, Time %20.14e, lv %d)\n", END_T, Time[lv], lv );
      Aux_Message( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n",
                   MPI_Rank, __FILE__, __LINE__, __FUNCTION__ );
      Aux_Message( stderr, "********************************************************************************\n" );
      MPI_Exit();
   }

   sprintf( dTime_Name[NdTime++], "%s", "End_Time" );


// 1.6 CRITERION SIX : user-defined criteria
// =============================================================================================================
   if ( OPT__DT_USER )
   {
      if ( Mis_GetTimeStep_User_Ptr != NULL )
      {
         dTime[NdTime] = dTime_dt * Mis_GetTimeStep_User_Ptr( lv, dTime_dt );
         sprintf( dTime_Name[NdTime++], "%s", "User" );
      }

      else
         Aux_Error( ERROR_INFO, "Mis_GetTimeStep_User_Ptr == NULL for OPT__DT_USER !!\n" );
   }


// 1.7 CRITERION SEVEN : phase rotation ##ELBDM ONLY##
// =============================================================================================================
#  if ( MODEL == ELBDM )
   if ( DT__PHASE != 0.0 )
   {
      dTime[NdTime] = dTime_dt * ELBDM_GetTimeStep_Phase( lv );
      sprintf( dTime_Name[NdTime++], "%s", "ELBDM_Phase" );
   }
#  endif


// 1.8 CRITERION EIGHT : particle evolution
// =============================================================================================================
#  ifdef PARTICLE
   Par_GetTimeStep_VelAcc( dTime[NdTime], dTime[NdTime+1], lv );

   dTime[NdTime] *= dTime_dt;
   sprintf( dTime_Name[NdTime++], "%s", "Par_Vel" );

   if ( DT__PARACC > 0.0 ) {
   dTime[NdTime] *= dTime_dt;
   sprintf( dTime_Name[NdTime++], "%s", "Par_Acc" ); }
#  endif



// 2. get the minimum time-step from all criteria
// =============================================================================================================
// 2.1 loop over all dt criteria
   double dTime_min = HUGE_NUMBER;

   for (int t=0; t<NdTime; t++)  dTime_min = fmin( dTime_min, dTime[t] );


// 2.2 record the physical time interval
// --> here we do not consider the time-step required to synchronize with the parent level (i.e., dTime_SyncFaLv)
   dTime_AllLv[lv] = dTime_min;


// 2.3 synchronize with the parent level
// --> increase dt at the current level by a small factor in order to synchronize with the parent level
   if ( OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE )
   {
      if ( lv > 0 )
      {
         if ( dTime_SyncFaLv <= 0.0 )
            Aux_Error( ERROR_INFO, "dTime_SyncFaLv (%20.14e) <= 0.0, something is wrong !!\n", dTime_SyncFaLv );

         if ( (1.0+DT__SYNC_PARENT_LV)*dTime_min >= dTime_SyncFaLv )    dTime_min = dTime_SyncFaLv;
      }

      dTime[NdTime] = dTime_SyncFaLv;
      sprintf( dTime_Name[NdTime++], "%s", "Sync_FaLv" );
   }


// 2.4 synchronize with the children level
// --> reduce dt at the current level by a small factor in order to synchronize with the children level easier
// --> assuming that dt at the children level is already known and won't change much in the next sub-step
// --> this could remove the additional sub-steps at the children level required to synchronize with this level
//     (these additional sub-steps usually have much smaller time-steps compared to the CFL condition)
// --> note that we do not apply this adjustment when this level is going to synchronize with its parent level
//     (i.e., dTime_min == dTime_SyncFaLv)
   if ( OPT__DT_LEVEL == DT_LEVEL_FLEXIBLE  &&  DT__SYNC_CHILDREN_LV > 0.0 )
   {
      const bool   Try2SyncSon     = ( lv < TOP_LEVEL  &&  NPatchTotal[lv+1] > 0  &&  dTime_min != dTime_SyncFaLv );
      const double dTime_SyncSonLv = ( Try2SyncSon ) ? 2.0*dTime_AllLv[lv+1] : NULL_REAL;

      if ( Try2SyncSon  &&  dTime_min > dTime_SyncSonLv  &&  dTime_min*(1.0-DT__SYNC_CHILDREN_LV) < dTime_SyncSonLv )
         dTime_min = dTime_SyncSonLv;

      dTime[NdTime] = dTime_SyncSonLv;
      sprintf( dTime_Name[NdTime++], "%s", "Sync_SonLv" );
   }


// 2.5 apply a maximum allowed dt
   if ( DT__MAX >= 0.0 )    dTime_min = MIN( dTime_min, DT__MAX );


// 2.6 reduce dt for AUTO_REDUCE_DT
// --> must do this AFTER checking all other dt criteria
   if ( AUTO_REDUCE_DT )   dTime_min *= AutoReduceDtCoeff;



// 3. record the dt info
// =============================================================================================================
   if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )
   {
      if ( FirstTime  &&  Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      FILE *File = fopen( FileName, "a" );

//    header
      if ( FirstTime )
      {
         fprintf( File, "#%3s  %8s  %8s  %13s  %13s  %13s", "Lv", "Step", "Counter", "TimeOld", "TimeNew", "dTime" );

#        ifdef COMOVING
         fprintf( File, "  %13s", "dTime_dt" );
#        endif

         for (int t=0; t<NdTime; t++)
         fprintf( File, "  %13s", dTime_Name[t] );

         if ( AUTO_REDUCE_DT )
         fprintf( File, "  %13s", "AutoRedDt" );

         fprintf( File, "\n" );
      }

//    dt info
      fprintf( File, "%4d  %8ld  %8ld  %13.7e  %13.7e  %13.7e",
               lv, Step, AdvanceCounter[lv], Time[lv], Time[lv]+dTime_min, dTime_min );

#     ifdef COMOVING
      fprintf( File, "  %13.7e", dTime_dt );
#     endif

      for (int t=0; t<NdTime; t++)
      fprintf( File, "  %13.7e", dTime[t] );

      if ( AUTO_REDUCE_DT )
      fprintf( File, "  %13.7e", AutoReduceDtCoeff );

      fprintf( File, "\n" );

      fclose( File );
   } // if ( OPT__RECORD_DT  &&  MPI_Rank == 0 )



// 4. verify time-step
// =============================================================================================================
   if ( dTime_min <= 0.0  ||  !Aux_IsFinite(dTime_min) )
      Aux_Error( ERROR_INFO, "incorrect time-step (dTime_min = %20.14e) !!\n", dTime_min );

// time synchronization may fail if dt/t is close to the round-off error limit
// --> e.g., temporal interpolation in Prepare_PatchData() may fail
   if (  Mis_CompareRealValue( Time[lv]+dTime_min, Time[lv], NULL, false )  )
      Aux_Error( ERROR_INFO, "time-step is too small (lv=%d, dt=%20.14e, t=%20.14e, dt/t=%20.14e, max_error=%20.14e) !!\n",
                 lv, dTime_min, Time[lv], dTime_min/Time[lv], MAX_ERROR_DBL );



   delete [] dTime_Name;
   delete [] dTime;
   FirstTime = false;

   return dTime_min;

} // FUNCTION : Mis_GetTimeStep
