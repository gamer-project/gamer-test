#include "GAMER.h"

#ifdef TIMING
extern Timer_t *Timer_dt         [NLEVEL];
extern Timer_t *Timer_Flu_Advance[NLEVEL];
extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_Src_Advance[NLEVEL];
extern Timer_t *Timer_Che_Advance[NLEVEL];
extern Timer_t *Timer_SF         [NLEVEL];
extern Timer_t *Timer_FB_Advance [NLEVEL];
extern Timer_t *Timer_FixUp      [NLEVEL];
extern Timer_t *Timer_Flag       [NLEVEL];
extern Timer_t *Timer_Refine     [NLEVEL];
extern Timer_t *Timer_Lv         [NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][9];
extern Timer_t *Timer_Par_Update [NLEVEL][3];
extern Timer_t *Timer_Par_2Sib   [NLEVEL];
extern Timer_t *Timer_Par_2Son   [NLEVEL];
#endif

bool AutoReduceDt_Continue;

extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime );




//-------------------------------------------------------------------------------------------------------
// Function    :  EvolveLevel
// Description :  Advance all physical attributes at the target level
//
// Note        :  1. Invoked by Main()
//
// Parameter   :  lv         : Target refinement level
//                dTime_FaLv : Physical time interval at the parent level
//-------------------------------------------------------------------------------------------------------
void EvolveLevel( const int lv, const double dTime_FaLv )
{

#  ifdef TIMING
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Lv[lv]->Start();
#  endif


#  ifdef GRAVITY
   const bool UsePot            = ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT );
#  endif
#  ifdef PARTICLE
   const bool StoreAcc_Yes      = true;
   const bool StoreAcc_No       = false;
   const bool UseStoredAcc_Yes  = true;
   const bool UseStoredAcc_No   = false;
   const bool TimingSendPar_Yes = true;
#  endif

   double dTime_SoFar, dTime_SubStep, dt_SubStep, TimeOld, TimeNew, AutoReduceDtCoeff;


// reset the workload weighting at each level to be recorded later
   if ( lv == 0 ) {
      for (int TLv=0; TLv<NLEVEL; TLv++)  amr->NUpdateLv[TLv] = 0; }


// sub-step loop
   dTime_SoFar           = 0.0;
   AutoReduceDtCoeff     = 1.0;
   AutoReduceDt_Continue = AUTO_REDUCE_DT;   // if AUTO_REDUCE_DT is on, we will perform auto-dt correction
                                             // at least once if the fluid solver fails

// note that we have ensured "Time[lv] == Time[lv-1]" to avoid the round-off errors
   while (  ( lv == 0 && amr->NUpdateLv[lv] == 0 )  ||
            ( lv >  0 && Time[lv] < Time[lv-1] )  )
   {
//    1. calculate the evolution time-step
// ===============================================================================================
#     ifdef TIMING
      if ( OPT__TIMING_BARRIER )    MPI_Barrier( MPI_COMM_WORLD );
      Timer_dt[lv]->Start();
#     endif

      switch ( OPT__DT_LEVEL )
      {
         case ( DT_LEVEL_SHARED ):
            if ( lv == 0 ) {
               dTime_SubStep = HUGE_NUMBER;

               for (int TLv=0; TLv<NLEVEL; TLv++)
               dTime_SubStep = fmin(  Mis_GetTimeStep( TLv, NULL_REAL, AutoReduceDtCoeff ), dTime_SubStep  );
            }
            else
               dTime_SubStep = dTime_FaLv;
         break;

         case ( DT_LEVEL_DIFF_BY_2 ):
            if ( lv == 0 ) {
               dTime_SubStep = HUGE_NUMBER;

               for (int TLv=0; TLv<NLEVEL; TLv++)
               dTime_SubStep = fmin(  Mis_GetTimeStep( TLv, NULL_REAL, AutoReduceDtCoeff )*(1<<TLv), dTime_SubStep  );
            }
            else
               dTime_SubStep = 0.5*dTime_FaLv;
         break;

         case ( DT_LEVEL_FLEXIBLE ):
               dTime_SubStep = Mis_GetTimeStep( lv, dTime_FaLv-dTime_SoFar, AutoReduceDtCoeff );
         break;

         default:
            Aux_Error( ERROR_INFO, "unsupported OPT__DT_LEVEL (%d) !!\n", OPT__DT_LEVEL );
      }

      dt_SubStep = Mis_dTime2dt( Time[lv], dTime_SubStep );
      TimeOld    = Time[lv];
      TimeNew    = Time[lv] + dTime_SubStep;

//    synchronize the time array to avoid the round-off errors
//    --> it's better to use dTime_SoFar+dTime_SubStep instead of Time[lv]+dTime_SubStep to check the synchronization
//        to reduce the round-off errors (in case Time[lv]>>dTime_SubStep)
      if (  lv > 0  &&  Mis_CompareRealValue( dTime_SoFar+dTime_SubStep, dTime_FaLv, NULL, false )  )
      TimeNew    = Time[lv-1];

      if ( lv == 0 )
      {
         dTime_Base = dTime_SubStep;

         if ( MPI_Rank == 0 )
            Aux_Message( stdout, "Time: %13.7e -> %13.7e,   Step: %7ld -> %7ld,   dt_base: %13.7e\n",
                         Time[0], Time[0]+dTime_Base, Step, Step+1, dTime_Base );
      }

#     ifdef TIMING
      Timer_dt[lv]->Stop();
#     endif
// ===============================================================================================


//    2. fluid solver
// ===============================================================================================
      const int SaveSg_Flu = 1 - amr->FluSg[lv];
#     ifdef MHD
      const int SaveSg_Mag = 1 - amr->MagSg[lv];
#     else
      const int SaveSg_Mag = NULL_INT;
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Flu_AdvanceDt, counter = %8ld ... ", lv, AdvanceCounter[lv] );

      if ( false ) {}
      /*
      if ( OPT__OVERLAP_MPI )
      {
//       enable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( true );
#        endif

//       advance patches needed to be sent
         TIMING_FUNC(   Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Mag, true, true ),
                        Timer_Flu_Advance[lv],   TIMER_ON   );

#        pragma omp parallel sections num_threads(2)
         {
#           pragma omp section
            {
//             transfer data simultaneously
#              ifdef GRAVITY
               if ( OPT__SELF_GRAVITY )
               TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT,   NULL_INT, DATA_GENERAL, _DENS,  _NONE, Rho_ParaBuf, USELB_YES ),
                              Timer_GetBuf[lv][0],   TIMER_ON   );
#              else
               TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, SaveSg_Mag, NULL_INT, DATA_GENERAL, _TOTAL, _MAG,  Flu_ParaBuf, USELB_YES ),
                              Timer_GetBuf[lv][2],   TIMER_ON   );
#              endif
            }

#           pragma omp section
            {
//             advance patches not needed to be sent
               TIMING_FUNC(   Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Mag, true, false ),
                              Timer_Flu_Advance[lv],   TIMER_ON   );
            }
         } // OpenMP parallel sections

//       disable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( false );
#        endif
      } // if ( OPT__OVERLAP_MPI )
      */

      else
      {
         int FluStatus_AllRank;

         TIMING_FUNC(   FluStatus_AllRank = Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Mag, false, false ),
                        Timer_Flu_Advance[lv],   TIMER_ON   );

//       do nothing if AUTO_REDUCE_DT is disabled
         if ( AUTO_REDUCE_DT )
         {
            if ( FluStatus_AllRank == GAMER_SUCCESS )
            {
//             restore the original parameters
               AutoReduceDtCoeff     = 1.0;
               AutoReduceDt_Continue = true;
            }

            else
            {
//             reduce the time-step coefficient if allowed
               if ( AutoReduceDtCoeff >= AUTO_REDUCE_DT_FACTOR_MIN )
               {
                  AutoReduceDtCoeff    *= AUTO_REDUCE_DT_FACTOR;
                  AutoReduceDt_Continue = true;

                  if ( MPI_Rank == 0 )
                  {
                     Aux_Message( stderr, "WARNING : fluid solver failed (Lv %2d, counter %8ld) --> ", lv, AdvanceCounter[lv] );
                     Aux_Message( stderr, "reduce dt by %13.7e\n", AutoReduceDtCoeff );
                  }
               }

//             if the time-step coefficient becomes smaller than the given threshold, restore the original time-step
//             and apply floor values in Flu_Close()
               else
               {
                  const double AutoReduceDtCoeff_Failed = AutoReduceDtCoeff;

                  AutoReduceDtCoeff     = 1.0;     // restore the original dt
                  AutoReduceDt_Continue = false;   // trigger density/energy floors in Flu_Close()

                  if ( MPI_Rank == 0 )
                  {
                     Aux_Message( stderr, "WARNING : AUTO_REDUCE_DT failed (Lv %2d, counter %8ld, dt-coeff=%13.7e < min=%13.7e) !!\n",
                                  lv, AdvanceCounter[lv], AutoReduceDtCoeff_Failed, AUTO_REDUCE_DT_FACTOR_MIN );
                     Aux_Message( stderr, "          --> Apply floor values with the original dt as the last resort ...\n" );
                  }
               } // if ( AutoReduceDtCoeff >= AUTO_REDUCE_DT_FACTOR_MIN ) ... else ...

//             restart the sub-step while loop
               continue;

            } // if ( FluStatus_AllRank == GAMER_SUCCESS ) ... else ...
         } // if ( AUTO_REDUCE_DT )

      } // if ( OPT__OVERLAP_MPI ) ... else ...

      amr->FluSg    [lv]             = SaveSg_Flu;
      amr->FluSgTime[lv][SaveSg_Flu] = TimeNew;
#     ifdef MHD
      amr->MagSg    [lv]             = SaveSg_Mag;
      amr->MagSgTime[lv][SaveSg_Mag] = TimeNew;
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================


//    3. update particles (prediction for KDK) and exchange particles
// ===============================================================================================
#     ifdef PARTICLE
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (predict) %5s... ", lv, "" );

#     ifdef STORE_PAR_ACC
      TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_PRED,
                                         (amr->Par->Integ == PAR_INTEG_EULER) ? StoreAcc_Yes    : StoreAcc_No,
                                         (amr->Par->Integ == PAR_INTEG_EULER) ? UseStoredAcc_No : UseStoredAcc_Yes ),
                     Timer_Par_Update[lv][0],   TIMER_ON   );
#     else
      TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_PRED,
                                         StoreAcc_No, UseStoredAcc_No ),
                     Timer_Par_Update[lv][0],   TIMER_ON   );
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Par_PassParticle2Sibling %9s... ", lv, "" );

      TIMING_FUNC(   Par_PassParticle2Sibling( lv, TimingSendPar_Yes ),
                     Timer_Par_2Sib[lv],   TIMER_ON   );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif
// ===============================================================================================


//    4. Poisson + gravity solver
// ===============================================================================================
#     ifdef GRAVITY
      const int SaveSg_Pot = 1 - amr->PotSg[lv];

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Gra_AdvanceDt, counter = %8ld ... ", lv, AdvanceCounter[lv] );

      if ( lv == 0 )
         Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot, UsePot, true, false, false, true );

      else // lv > 0
      {
         if ( false ) {}
         /*
         if ( OPT__OVERLAP_MPI )
         {
//          enable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( true );
#           endif

//          advance patches needed to be sent
            TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                           UsePot, true, true, true, true ),
                           Timer_Gra_Advance[lv],   TIMER_ON   );

#           pragma omp parallel sections num_threads(2)
            {
#              pragma omp section
               {
//                transfer data simultaneously
                  if ( UsePot )
                  TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON,
                                                    _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][1],   TIMER_ON   );

                  TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, SaveSg_Mag, NULL_INT, DATA_GENERAL,
                                                    _TOTAL, _MAG,  Flu_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][2],   TIMER_ON   );
               }

#              pragma omp section
               {
//                advance patches not needed to be sent
                  TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                                                UsePot, true, true, false, true ),
                                 Timer_Gra_Advance[lv],   TIMER_ON   );
               }
            } // OpenMP parallel sections

//          disable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( false );
#           endif
         } // if ( OPT__OVERLAP_MPI )
         */

         else
         {
//          exchange the updated density field in the buffer patches for the Poisson solver
            if ( OPT__SELF_GRAVITY )
            TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, NULL_INT, DATA_GENERAL,
                                              _DENS, _NONE, Rho_ParaBuf, USELB_YES ),
                           Timer_GetBuf[lv][0],   TIMER_ON   );

            TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                                          UsePot, true, false, false, true ),
                           Timer_Gra_Advance[lv],   TIMER_ON   );

//          exchange the updated potential in the buffer patches
//          --> we will do this after all other operations (e.g., star formation) if OPT__MINIMIZE_MPI_BARRIER is adopted
//              --> assuming that all remaining operations do not need to access the potential in the buffer patches
//              --> one must enable both STORE_POT_GHOST and PAR_IMPROVE_ACC for this purpose
            if ( UsePot  &&  !OPT__MINIMIZE_MPI_BARRIER )
            TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON,
                                              _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                           Timer_GetBuf[lv][1],   TIMER_ON   );
         } // if ( OPT__OVERLAP_MPI ) ... else ...

         if ( UsePot )
         {
            amr->PotSg    [lv]             = SaveSg_Pot;
            amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;
         }

      } // if ( lv == 0 ) ... else ...

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif // #ifdef GRAVITY
// ===============================================================================================


//    5. correct particles velocity and send particles to lv+1
// ===============================================================================================
#     ifdef PARTICLE
      if ( amr->Par->Integ == PAR_INTEG_KDK )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (correct Lv %2d)... ", lv, lv );

#        ifdef STORE_PAR_ACC
         TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_CORR, StoreAcc_Yes, UseStoredAcc_No ),
                        Timer_Par_Update[lv][1],   TIMER_ON   );
#        else
         TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_CORR, StoreAcc_No,  UseStoredAcc_No ),
                        Timer_Par_Update[lv][1],   TIMER_ON   );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

         if ( lv > 0 )
         {
            if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
               Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (correct Lv %2d)... ", lv, lv-1 );

//          apply velocity correction for particles just travelling from lv to lv-1
#           ifdef STORE_PAR_ACC
            TIMING_FUNC(   Par_UpdateParticle( lv-1, TimeNew, TimeOld, PAR_UPSTEP_CORR, StoreAcc_Yes, UseStoredAcc_No ),
                           Timer_Par_Update[lv][2],   TIMER_ON   );
#           else
            TIMING_FUNC(   Par_UpdateParticle( lv-1, TimeNew, TimeOld, PAR_UPSTEP_CORR, StoreAcc_No,  UseStoredAcc_No ),
                           Timer_Par_Update[lv][2],   TIMER_ON   );
#           endif

            if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
         }
      }

//    pass particles to the children patches
//    --> we will do this later (just before the star-formation routines) if OPT__MINIMIZE_MPI_BARRIER is adopted
      if ( !OPT__MINIMIZE_MPI_BARRIER )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Par_PassParticle2Son %12s ... ", lv, "" );

         TIMING_FUNC(   Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_Yes, NULL_INT, NULL ),
                        Timer_Par_2Son[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      }
#     endif // #ifdef PARTICLE
// ===============================================================================================


//    6. additional physics
// ===============================================================================================

// *********************************
//    6-1. local source terms
// *********************************
      const int SaveSg_SrcFlu = SaveSg_Flu;  // save in the same Flu/MagSg
      const int SaveSg_SrcMag = SaveSg_Mag;

      if ( SrcTerms.Any )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Src_AdvanceDt, counter = %8ld ... ", lv, AdvanceCounter[lv] );

//###REVISE: we have assumed that Src_AdvanceDt() requires no ghost zones
         TIMING_FUNC(   Src_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_SrcFlu, SaveSg_SrcMag, false, false ),
                        Timer_Src_Advance[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      }


// *********************************
//    6-2. Grackle cooling/heating
// *********************************
#     ifdef SUPPORT_GRACKLE
      if ( GRACKLE_ACTIVATE )
      {
         const int SaveSg_Che = SaveSg_Flu;  // save in the same FluSg

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Grackle_AdvanceDt, counter = %4ld ... ", lv, AdvanceCounter[lv] );

//###REVISE: we have assumed that Grackle_AdvanceDt() requires no ghost zones
         TIMING_FUNC(   Grackle_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Che, false, false ),
                        Timer_Che_Advance[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // if ( GRACKLE_ACTIVATE )
#     endif // #ifdef SUPPORT_GRACKLE


// *********************************
//    6-3. star formation
// *********************************
#     ifdef PARTICLE
//    pass particles to the children patches here if OPT__MINIMIZE_MPI_BARRIER is adopted
//    --> do this before any star-formation and feedback routines so that particles always live in the leaf patches
      if ( OPT__MINIMIZE_MPI_BARRIER )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Par_PassParticle2Son %12s ... ", lv, "" );

         TIMING_FUNC(   Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_Yes, NULL_INT, NULL ),
                        Timer_Par_2Son[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      }
#     endif // #ifdef PARTICLE


#     ifdef STAR_FORMATION
      if ( SF_CREATE_STAR_SCHEME != SF_CREATE_STAR_SCHEME_NONE )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: SF_CreateStar, counter = %8ld ... ", lv, AdvanceCounter[lv] );

//###REVISE: we have assumed that SF_CreateStar() requires no ghost zones
         TIMING_FUNC(   SF_CreateStar( lv, TimeNew, dt_SubStep ),
                        Timer_SF[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // if ( SF_CREATE_STAR_SCHEME != SF_CREATE_STAR_SCHEME_NONE )
#     endif // #ifdef STAR_FORMATION


// *********************************
//    6-4. feedback
// *********************************
#     ifdef FEEDBACK
      const int SaveSg_FBFlu = SaveSg_Flu;   // save in the same Flu/MagSg
      const int SaveSg_FBMag = SaveSg_Mag;

      if ( true )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: FB_AdvanceDt, counter = %9ld ... ", lv, AdvanceCounter[lv] );

//###REVISE: we have assumed that FB_AdvanceDt() requires no ghost zones
         TIMING_FUNC(   FB_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_FBFlu, SaveSg_FBMag ),
                        Timer_FB_Advance[lv],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      }
#     endif // #ifdef FEEDBACK

// ===============================================================================================


//    7. overwrite the fluid fields
// ===============================================================================================
      if ( OPT__RESET_FLUID )
      {
//       use the same timer as the fluid solver for now
         if ( Flu_ResetByUser_API_Ptr != NULL )
         {
            TIMING_FUNC(   Flu_ResetByUser_API_Ptr( lv, SaveSg_Flu, TimeNew ),
                           Timer_Flu_Advance[lv],   TIMER_ON   );
         }

         else
            Aux_Error( ERROR_INFO, "Flu_ResetByUser_API_Ptr == NULL for OPT__RESET_FLUID !!\n" );
      }
// ===============================================================================================


//    exchange the updated fluid field in the buffer patches
      TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, SaveSg_Mag, NULL_INT, DATA_GENERAL,
                                        _TOTAL, _MAG, Flu_ParaBuf, USELB_YES ),
                     Timer_GetBuf[lv][2],   TIMER_ON   );

//    exchange the updated potential in the buffer patches here if OPT__MINIMIZE_MPI_BARRIER is adopted
#     ifdef GRAVITY
      if ( lv > 0  &&  UsePot  &&  OPT__MINIMIZE_MPI_BARRIER )
      TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON,
                                        _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                     Timer_GetBuf[lv][1],   TIMER_ON   );
#     endif


      dTime_SoFar       += dTime_SubStep;
      Time_Prev     [lv] = TimeOld;
      Time          [lv] = TimeNew;
      AdvanceCounter[lv] ++;
      amr->NUpdateLv[lv] ++;

      if ( AdvanceCounter[lv] >= __LONG_MAX__ )    Aux_Message( stderr, "WARNING : AdvanceCounter overflow !!\n" );


      if ( lv != TOP_LEVEL  &&  NPatchTotal[lv+1] != 0 )
      {

//       7. enter the next refinement level
// ===============================================================================================
#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Lv[lv]->Stop();
#        endif

         EvolveLevel( lv+1, dTime_SubStep );

#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Lv[lv]->Start();
#        endif
// ===============================================================================================


//       8. correct the data at the current level with the data at the next finer level
// ===============================================================================================
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flu_FixUp %24s... ", lv, "" );

//       8-1. use the average data on fine grids to correct the coarse-grid data
         if ( OPT__FIXUP_RESTRICT )
         {
            TIMING_FUNC(   Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv],
                                               NULL_INT, NULL_INT, _TOTAL, _MAG ),
                           Timer_FixUp[lv],   TIMER_ON   );

#           ifdef LOAD_BALANCE
            TIMING_FUNC(   LB_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_RESTRICT,
                                             _TOTAL, _MAG, NULL_INT ),
                           Timer_GetBuf[lv][7],   TIMER_ON   );
#           endif
         }

//       8-2. use the fine-grid electric field on the coarse-fine boundaries to correct the coarse-grid magnetic field
#        ifdef MHD
         if ( OPT__FIXUP_ELECTRIC )
         {
#           ifdef LOAD_BALANCE
            TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, NULL_INT, COARSE_FINE_ELECTRIC,
                                              _NONE, _NONE, NULL_INT, USELB_YES ),
                           Timer_GetBuf[lv][6],   TIMER_ON   );
#           endif

            TIMING_FUNC(   MHD_FixUp_Electric( lv ),
                           Timer_FixUp[lv],   TIMER_ON   );
         }
#        endif

//       8-3. use the fine-grid fluxes across the coarse-fine boundaries to correct the coarse-grid data
//            --> apply AFTER other fix-up operations since it will check negative pressure as well
//                (which requires the coarse-grid B field updated by Flu_FixUp_Restrict() and MHD_FixUp_Electric())
         if ( OPT__FIXUP_FLUX )
         {
#           ifdef LOAD_BALANCE
            TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, NULL_INT, COARSE_FINE_FLUX,
                                              _FLUX_TOTAL, _NONE, NULL_INT, USELB_YES ),
                           Timer_GetBuf[lv][6],   TIMER_ON   );
#           endif

            TIMING_FUNC(   Flu_FixUp_Flux( lv ),
                           Timer_FixUp[lv],   TIMER_ON   );
         }

//       8-4. exchange the updated data
#        ifdef MHD
         if ( OPT__FIXUP_FLUX  ||  OPT__FIXUP_RESTRICT  ||  OPT__FIXUP_ELECTRIC )
#        else
         if ( OPT__FIXUP_FLUX  ||  OPT__FIXUP_RESTRICT )
#        endif
         TIMING_FUNC(   Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_AFTER_FIXUP,
                                           _TOTAL, _MAG, Flu_ParaBuf, USELB_YES  ),
                        Timer_GetBuf[lv][3],   TIMER_ON   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================

      } // if ( lv != TOP_LEVEL  &&  NPatchTotal[lv+1] != 0 )


//    9. flag the current level and create patches at the next finer level
// ===============================================================================================
      if ( lv != TOP_LEVEL  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
      {
//       flag
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flag %29s... ", lv, "" );

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Flag_Real( lv, USELB_YES ),       Timer_Flag[lv],   TIMER_ON   );
#        else
         TIMING_FUNC(   Flag_Real( lv, USELB_NO ),        Timer_Flag[lv],   TIMER_ON   );

         TIMING_FUNC(   MPI_ExchangeBoundaryFlag( lv ),   Timer_Flag[lv],   TIMER_ON   );

         TIMING_FUNC(   Flag_Buffer( lv ),                Timer_Flag[lv],   TIMER_ON   );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


//       refine
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Refine %27s... ", lv, "" );

         TIMING_FUNC(   Refine( lv, USELB_YES ),
                        Timer_Refine[lv],   TIMER_ON   );

         Time          [lv+1]                     = Time[lv];
         amr->FluSgTime[lv+1][ amr->FluSg[lv+1] ] = Time[lv];
#        ifdef MHD
         amr->MagSgTime[lv+1][ amr->MagSg[lv+1] ] = Time[lv];
#        endif
#        ifdef GRAVITY
         if ( UsePot )
         amr->PotSgTime[lv+1][ amr->PotSg[lv+1] ] = Time[lv];
#        endif

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_AFTER_REFINE,
                                           _TOTAL, _MAG, Flu_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][4],   TIMER_ON   );
#        ifdef GRAVITY
         if ( UsePot )
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_AFTER_REFINE,
                                           _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][5],   TIMER_ON   );
#        endif
#        endif // #ifdef LOAD_BALANCE

         TIMING_FUNC(   Buf_GetBufferData( lv+1, amr->FluSg[lv+1], amr->MagSg[lv+1], NULL_INT, DATA_AFTER_REFINE,
                                           _TOTAL, _MAG, Flu_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][4],   TIMER_ON   );
#        ifdef GRAVITY
         if ( UsePot )
         TIMING_FUNC(   Buf_GetBufferData( lv+1, NULL_INT, NULL_INT, amr->PotSg[lv+1], POT_AFTER_REFINE,
                                           _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][5],   TIMER_ON   );
#        endif

//       must call Poi_StorePotWithGhostZone() AFTER collecting potential for buffer patches
#        ifdef STORE_POT_GHOST
         if ( UsePot )
         TIMING_FUNC(   Poi_StorePotWithGhostZone( lv+1, amr->PotSg[lv+1], false ),
                        Timer_Refine[lv],   TIMER_ON   );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

         if ( OPT__PATCH_COUNT == 2 )  Aux_Record_PatchCount();
#        ifdef PARTICLE
         if ( OPT__PARTICLE_COUNT == 2 )  Par_Aux_Record_ParticleCount();
#        endif

      } // if ( lv != TOP_LEVEL  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
// ===============================================================================================

   } // while()


#  ifdef TIMING
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Lv[lv]->Stop();
#  endif

} // FUNCTION : EvolveLevel

