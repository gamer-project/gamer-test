#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp_Restrict
// Description :  Replace the data at level "FaLv" by the average data at level "FaLv+1"
//
// Note        :  1. Use the input parameters "TVarCC" and "TVarFC" to determine the target cell- and
//                   face-centered variables, respectively
//                2. For MHD, this routine currently always restrict all three B field components
//                   --> Do not distinguish _MAGX, _MAGY, _MAGZ, and _MAG in TVarFC
//                3. Invoked by EvolveLevel()
//
// Parameter   :  FaLv     : Target refinement level at which the data are going to be replaced
//                SonFluSg : Fluid sandglass at level "FaLv+1"
//                FaFluSg  : Fluid sandglass at level "FaLv"
//                SonMagSg : B field sandglass at level "FaLv+1"
//                FaMagSg  : B field sandglass at level "FaLv"
//                SonPotSg : Potential sandglass at level "FaLv+1"
//                FaPotSg  : Potential sandglass at level "FaLv"
//                TVarCC   : Target cell-centered variables
//                           --> Supported variables in different models:
//                               HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY,[, _POTE]
//                               ELBDM : _DENS, _REAL, _IMAG, [, _POTE]
//                           --> _FLUID, _PASSIVE, and _TOTAL apply to all models
//                TVarFC   : Target face-centered variables
//                            --> Supported variables in different models:
//                                HYDRO with MHD : _MAGX, _MAGY, _MAGZ, _MAG
//                                ELBDM          : none
//                            --> But it currently does not distinguish _MAGX, _MAGY, _MAGZ, and _MAG
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonMagSg, const int FaMagSg,
                         const int SonPotSg, const int FaPotSg, const long TVarCC, const long TVarFC )
{

   const int SonLv = FaLv + 1;

// check
   if ( FaLv < 0  ||  FaLv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );

   if ( FaLv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : applying \"%s\" to the maximum level is meaningless !! \n", __FUNCTION__ );
      return;
   }

   if (  ( TVarCC & _TOTAL )  &&  ( SonFluSg != 0 && SonFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonFluSg", SonFluSg );

   if (  ( TVarCC & _TOTAL )  &&  ( FaFluSg != 0 && FaFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaFluSg", FaFluSg );

#  ifdef GRAVITY
   if (  ( TVarCC & _POTE )  &&  ( SonPotSg != 0 && SonPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonPotSg", SonPotSg );

   if (  ( TVarCC & _POTE )  &&  ( FaPotSg != 0 && FaPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPotSg", FaPotSg );
#  endif

#  ifdef MHD
   if (  ( TVarFC & _MAG )  &&  TVarFC != _MAG  )
      Aux_Error( ERROR_INFO, "must work on all three magnetic components at once !!\n" );

   if (  ( TVarFC & _MAG )  &&  ( SonMagSg != 0 && SonMagSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonMagSg", SonMagSg );

   if (  ( TVarFC & _MAG )  &&  ( FaMagSg != 0 && FaMagSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaMagSg", FaMagSg );
#  endif


// nothing to do if there are no real patches at lv+1
   if ( amr->NPatchComma[SonLv][1] == 0 )    return;

// check the synchronization
   Mis_CompareRealValue( Time[FaLv], Time[SonLv], __FUNCTION__, true );


   const bool ResFlu  = TVarCC & _TOTAL;
#  ifdef GRAVITY
   const bool ResPot  = TVarCC & _POTE;
#  else
   const bool ResPot  = false;
#  endif
#  ifdef MHD
   const bool ResMag  = TVarFC & _MAG;
#  else
   const bool ResMag  = false;
#  endif
   const int PS1_half = PS1 / 2;


// determine the components to be restricted (TFluVarIdx : target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
   int NFluVar=0, TFluVarIdxList[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVarCC & (1L<<v) )    TFluVarIdxList[ NFluVar ++ ] = v;


// return immediately if no target variable is found
   if ( NFluVar == 0  &&  !ResPot  &&  !ResMag )
   {
      Aux_Message( stderr, "WARNING : no target variable is found !!\n" );
      return;
   }


// restrict
#  pragma omp parallel for schedule( runtime )
   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)
   {
      const int FaPID = amr->patch[0][SonLv][SonPID0]->father;

//    check
#     ifdef GAMER_DEBUG
      if ( FaPID < 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch (FaPID = %d) !!\n",
                    SonLv, SonPID0, FaPID );

      if ( ResFlu  &&  amr->patch[FaFluSg][FaLv][FaPID]->fluid == NULL )
         Aux_Error( ERROR_INFO, "FaFluSg %d, FaLv %d, FaPID %d has no fluid array allocated !!\n",
                    FaFluSg, FaLv, FaPID );

#     ifdef GRAVITY
      if ( ResPot  &&  amr->patch[FaPotSg][FaLv][FaPID]->pot == NULL )
         Aux_Error( ERROR_INFO, "FaPotSg %d, FaLv %d, FaPID %d has no potential array allocated !!\n",
                    FaPotSg, FaLv, FaPID );
#     endif

#     ifdef MHD
      if ( ResMag  &&  amr->patch[FaMagSg][FaLv][FaPID]->magnetic == NULL )
         Aux_Error( ERROR_INFO, "FaMagSg %d, FaLv %d, FaPID %d has no B field array allocated !!\n",
                    FaMagSg, FaLv, FaPID );
#     endif
#     endif // #ifdef GAMER_DEBUG


//    loop over eight sons
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int SonPID = SonPID0 + LocalID;
         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1_half );
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1_half );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1_half );

//       check
#        ifdef GAMER_DEBUG
         if ( ResFlu  &&  amr->patch[SonFluSg][SonLv][SonPID]->fluid == NULL )
            Aux_Error( ERROR_INFO, "SonFluSg %d, SonLv %d, SonPID %d has no fluid array allocated !!\n",
                       SonFluSg, SonLv, SonPID );

#        ifdef GRAVITY
         if ( ResPot  &&  amr->patch[SonPotSg][SonLv][SonPID]->pot == NULL )
            Aux_Error( ERROR_INFO, "SonPotSg %d, SonLv %d, SonPID %d has no potential array allocated !!\n",
                       SonPotSg, SonLv, SonPID );
#        endif

#        ifdef MHD
         if ( ResMag  &&  amr->patch[SonMagSg][SonLv][SonPID]->magnetic == NULL )
            Aux_Error( ERROR_INFO, "SonMagSg %d, SonLv %d, SonPID %d has no B field array allocated !!\n",
                       SonMagSg, SonLv, SonPID );
#        endif
#        endif // #ifdef GAMER_DEBUG


//       restrict the fluid data
         if ( ResFlu )
         for (int v=0; v<NFluVar; v++)
         {
            const int TFluVarIdx = TFluVarIdxList[v];
            const real (*SonPtr)[PS1][PS1] = amr->patch[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx];
                  real (* FaPtr)[PS1][PS1] = amr->patch[ FaFluSg][ FaLv][ FaPID]->fluid[TFluVarIdx];

            int ii, jj, kk, I, J, K, Ip, Jp, Kp;

            for (int k=0; k<PS1_half; k++)  {  K = k*2;  Kp = K+1;  kk = k + Disp_k;
            for (int j=0; j<PS1_half; j++)  {  J = j*2;  Jp = J+1;  jj = j + Disp_j;
            for (int i=0; i<PS1_half; i++)  {  I = i*2;  Ip = I+1;  ii = i + Disp_i;

               FaPtr[kk][jj][ii] = 0.125*( SonPtr[K ][J ][I ] + SonPtr[K ][J ][Ip] +
                                           SonPtr[K ][Jp][I ] + SonPtr[Kp][J ][I ] +
                                           SonPtr[K ][Jp][Ip] + SonPtr[Kp][Jp][I ] +
                                           SonPtr[Kp][J ][Ip] + SonPtr[Kp][Jp][Ip] );
            }}}
         } // if ( ResFlu )


//       restrict the potential data
#        ifdef GRAVITY
         if ( ResPot )
         {
            const real (*SonPtr)[PS1][PS1] = amr->patch[SonPotSg][SonLv][SonPID]->pot;
                  real (* FaPtr)[PS1][PS1] = amr->patch[ FaPotSg][ FaLv][ FaPID]->pot;

            int ii, jj, kk, I, J, K, Ip, Jp, Kp;

            for (int k=0; k<PS1_half; k++)  {  K = k*2;  Kp = K+1;  kk = k + Disp_k;
            for (int j=0; j<PS1_half; j++)  {  J = j*2;  Jp = J+1;  jj = j + Disp_j;
            for (int i=0; i<PS1_half; i++)  {  I = i*2;  Ip = I+1;  ii = i + Disp_i;

               FaPtr[kk][jj][ii] = 0.125*( SonPtr[K ][J ][I ] + SonPtr[K ][J ][Ip] +
                                           SonPtr[K ][Jp][I ] + SonPtr[Kp][J ][I ] +
                                           SonPtr[K ][Jp][Ip] + SonPtr[Kp][Jp][I ] +
                                           SonPtr[Kp][J ][Ip] + SonPtr[Kp][Jp][Ip] );
            }}}
         } // if ( ResPot )
#        endif // #ifdef GRAVITY


//       restrict the magnetic field
//       --> currently it always works on all three B field components
//###OPTIMIZATION: coarse-grid B field on the son patch boundaries (within the same patch group)
//                 is restricted twice by two adjacent son patches
#        ifdef MHD
         if ( ResMag )
         {
            int idx_fa, idx_son0, I, J, K;

//          Bx
            const real *SonBx = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[0];
                  real * FaBx = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[0];

            for (int k=0; k<PS1_half;   k++)  {  K = k*2;
            for (int j=0; j<PS1_half;   j++)  {  J = j*2;
            for (int i=0; i<PS1_half+1; i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1P1, PS1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1P1, PS1 );

               FaBx[idx_fa] = 0.25*( SonBx[ idx_son0                     ] +
                                     SonBx[ idx_son0 + PS1P1             ] +
                                     SonBx[ idx_son0 + PS1P1*PS1         ] +
                                     SonBx[ idx_son0 + PS1P1*PS1 + PS1P1 ] );
            }}}

//          By
            const real *SonBy = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[1];
                  real * FaBy = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[1];

            for (int k=0; k<PS1_half;   k++)  {  K = k*2;
            for (int j=0; j<PS1_half+1; j++)  {  J = j*2;
            for (int i=0; i<PS1_half;   i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1, PS1P1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1, PS1P1 );

               FaBy[idx_fa] = 0.25*( SonBy[ idx_son0                 ] +
                                     SonBy[ idx_son0 + 1             ] +
                                     SonBy[ idx_son0 + PS1P1*PS1     ] +
                                     SonBy[ idx_son0 + PS1P1*PS1 + 1 ] );
            }}}

//          Bz
            const real *SonBz = amr->patch[SonMagSg][SonLv][SonPID]->magnetic[2];
                  real * FaBz = amr->patch[ FaMagSg][ FaLv][ FaPID]->magnetic[2];

            for (int k=0; k<PS1_half+1; k++)  {  K = k*2;
            for (int j=0; j<PS1_half;   j++)  {  J = j*2;
            for (int i=0; i<PS1_half;   i++)  {  I = i*2;

               const int idx_fa   = IDX321( i+Disp_i, j+Disp_j, k+Disp_k, PS1, PS1 );
               const int idx_son0 = IDX321( I,        J,        K,        PS1, PS1 );

               FaBz[idx_fa] = 0.25*( SonBz[ idx_son0           ] +
                                     SonBz[ idx_son0 + 1       ] +
                                     SonBz[ idx_son0 + PS1     ] +
                                     SonBz[ idx_son0 + PS1 + 1 ] );
            }}}
         } // if ( ResMag )
#        endif // ifdef MHD
      } // for (int LocalID=0; LocalID<8; LocalID++)


//    apply the same B field restriction to the data of father-sibling patches on the coarse-fine boundaries
#     ifdef MHD
      for (int s=0; s<6; s++)
      {
         const int FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#        ifdef GAMER_DEBUG
         if ( FaSibPID == -1 )
            Aux_Error( ERROR_INFO, "FaSibPID == -1 (FaLv %d, FaPID %d, s %d, SonPID0 %d) !!\n", FaLv, FaPID, s, SonPID0 );
#        endif

//       skip father patches adjacent to non-periodic boundaries
         if ( FaSibPID < -1 )    continue;

//       find the coarse-fine boundaries and copy data
         if ( amr->patch[0][FaLv][FaSibPID]->son == -1 )    MHD_CopyPatchInterfaceBField( FaLv, FaPID, s, FaMagSg );
      }
#     endif // #ifdef MHD


//    check the minimum pressure/internal energy and, when the dual-energy formalism is adopted, ensure the consistency between
//    pressure, total energy density, and the dual-energy variable
#     if ( MODEL == HYDRO && !defined SRHD )
//    apply this correction only when preparing all fluid variables or magnetic field
#     ifdef MHD
      if (  ( TVarCC & _TOTAL ) == _TOTAL  ||  ResMag  )
#     else
      if (  ( TVarCC & _TOTAL ) == _TOTAL  )
#     endif
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
//       compute magnetic energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( FaLv, FaPID, i, j, k, FaMagSg );
#        else
         const real Emag = NULL_REAL;
#        endif

#        ifdef DUAL_ENERGY
//       here we ALWAYS use the dual-energy variable to correct the total energy density
//       --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//           the runtime parameter DUAL_ENERGY_SWITCH here
         const bool CheckMinPres_Yes = true;
         const real UseEnpy2FixEngy  = HUGE_NUMBER;
         char dummy;    // we do not record the dual-energy status here

         Hydro_DualEnergyFix( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                              amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENPY][k][j][i],
                              dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2], CheckMinPres_Yes, MIN_PRES,
                              UseEnpy2FixEngy, Emag );

#        else // #ifdef DUAL_ENERGY

//       actually it might not be necessary to check the minimum internal energy here
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i]
            = Hydro_CheckMinEintInEngy( amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMX][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMY][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[MOMZ][k][j][i],
                                        amr->patch[FaFluSg][FaLv][FaPID]->fluid[ENGY][k][j][i],
                                        MIN_EINT, Emag );
#        endif // #ifdef DUAL_ENERGY ... else ...
      } // i,j,k
#     endif // #if ( MODEL == HYDRO )


//    rescale real and imaginary parts to get the correct density in ELBDM
#     if ( MODEL == ELBDM )
      real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

      if (  ( TVarCC & _DENS )  &&  ( TVarCC & _REAL )  &&  (TVarCC & _IMAG )  )
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         Real      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i];
         Imag      = amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
         {
            amr->patch[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i] = (real)0.0;
            Rescale = (real)0.0;
         }
         else
            Rescale = SQRT( Rho_Corr/Rho_Wrong );

         amr->patch[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i] *= Rescale;
         amr->patch[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i] *= Rescale;
      }
#     endif

   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)

} // FUNCTION : Flu_FixUp_Restrict
