#include "GAMER.h"
#include "CUFLU.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp_Flux
// Description :  Use the fine-grid fluxes across coarse-fine boundaries to correct the coarse-grid data
//
// Note        :  1. Boundary fluxes from the neighboring ranks must be received in advance by invoking
//                   Buf_GetBufferData()
//                2. Invoked by EvolveLevel()
//
// Parameter   :  lv : Target coarse level
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp_Flux( const int lv )
{

   const real Const[6]   = { real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]),
                             real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]),
                             real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]) };
   const int  FluSg      = amr->FluSg[lv];
#  ifdef MHD
   const int  MagSg      = amr->MagSg[lv];
#  endif
   const int  Offset[6]  = { 0, PS1-1, 0, (PS1-1)*PS1, 0, (PS1-1)*SQR(PS1) }; // x=0/PS1-1, y=0/PS1-1, z=0/PS1-1 faces
   const int  didx[3][2] = { PS1, SQR(PS1), 1, SQR(PS1), 1, PS1 };

//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
// when enabling cooling, we want to fix the specific internal energy so that it won't be corrected by the flux fix-up operation
// --> it doesn't make much sense to correct the total energy density, pressure, or dual-energy variable when cooling is
//     adopted since the total energy is not conserved anymore
// --> experiments show that fixing the specific internal energy instead of pressure works better
//     --> fixing pressure is found to lead to extremely small dt
//     --> fixing specific internal energy works better since it is directly proportional to the sound speed square
#  ifdef SUPPORT_GRACKLE
   const bool FixSEint   = GRACKLE_ACTIVATE;
#  else
   const bool FixSEint   = false;
#  endif
*/


// check
#  ifdef GAMER_DEBUG

   if ( !amr->WithFlux )
      Aux_Error( ERROR_INFO, "amr->WithFlux is off -> no flux array is allocated for OPT__FIXUP_FLUX !!\n" );

#  if ( MODEL == ELBDM )

#  ifndef CONSERVE_MASS
   Aux_Error( ERROR_INFO, "CONSERVE_MASS is not turned on in the Makefile for the option OPT__FIXUP_FLUX !!\n" );
#  endif

#  if ( NFLUX_TOTAL != 1 )
   Aux_Error( ERROR_INFO, "NFLUX_TOTAL (%d) != 1 for the option OPT__FIXUP_FLUX !!\n", NFLUX_TOTAL );
#  endif

#  if ( DENS != 0 )
   Aux_Error( ERROR_INFO, "DENS (%d) != 0 for the option OPT__FIXUP_FLUX !!\n", DENS );
#  endif

#  if ( FLUX_DENS != 0 )
   Aux_Error( ERROR_INFO, "FLUX_DENS (%d) != 0 for the option OPT__FIXUP_FLUX !!\n", FLUX_DENS );
#  endif


// if "NCOMP_TOTAL != NFLUX_TOTAL", one must specify how to correct cell data from the flux arrays
// --> specifically, how to map different flux variables to fluid active/passive variables
// --> for ELBDM, we have assumed that FLUX_DENS == DENS == 0 and NFLUX_TOTAL == 1
#  elif ( NCOMP_TOTAL != NFLUX_TOTAL )
#     error : NCOMP_TOTAL != NFLUX_TOTAL (one must specify how to map flux variables to fluid active/passive variables) !!

#  endif // #if ( MODEL == ELBDM ) ... elif ...

#  endif // #ifdef GAMER_DEBUG


#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    1. sum up the coarse-grid and fine-grid fluxes for bitwise reproducibility
#     ifdef BIT_REP_FLUX
      for (int s=0; s<6; s++)
      {
         real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux[s];

         if ( FluxPtr != NULL )
         {
            for (int v=0; v<NFLUX_TOTAL; v++)
            for (int m=0; m<PS1; m++)
            for (int n=0; n<PS1; n++)
               FluxPtr[v][m][n] += amr->patch[0][lv][PID]->flux_bitrep[s][v][m][n];
         }
      }
#     endif


//    2. correct fluid variables by the difference between the coarse-grid and fine-grid fluxes
//    loop over all six faces of a given patch
      for (int s=0; s<6; s++)
      {
//       skip the faces not adjacent to the coarse-fine boundaries
         const real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux[s];
         if ( FluxPtr == NULL  )  continue;


//       set the pointers to the target face
         real *FluidPtr1D0[NCOMP_TOTAL], *FluidPtr1D[NCOMP_TOTAL];
         for (int v=0; v<NCOMP_TOTAL; v++)   FluidPtr1D0[v] = amr->patch[FluSg][lv][PID]->fluid[v][0][0] + Offset[s];
#        ifdef DUAL_ENERGY
         const char *DE_StatusPtr1D0 = amr->patch[0][lv][PID]->de_status[0][0] + Offset[s];
#        endif


//       set the array index strides
         const int didx_m = didx[s/2][1];
         const int didx_n = didx[s/2][0];


//       loop over all cells on a given face
         for (int m=0; m<PS1; m++)
         {
            for (int v=0; v<NCOMP_TOTAL; v++)   FluidPtr1D[v] = FluidPtr1D0[v] + m*didx_m;
#           ifdef DUAL_ENERGY
            const char *DE_StatusPtr1D = DE_StatusPtr1D0 + m*didx_m;
#           endif

            for (int n=0; n<PS1; n++)
            {
//             from now on we also correct cells updated by either the minimum pressure threshold or the 1st-order-flux correction
//             --> note that we have stored the 1st-order fluxes across the coarse-fine boundaries in Flu_Close()
               /*
#              ifdef DUAL_ENERGY
               if ( *DE_StatusPtr1D == DE_UPDATED_BY_MIN_PRES  ||  *DE_StatusPtr1D == DE_UPDATED_BY_1ST_FLUX )    continue;
#              endif
               */


//             calculate the corrected results
//             --> do NOT **store** these results yet since we want to skip the cells with unphysical results
               real CorrVal[NFLUX_TOTAL];    // values after applying the flux correction
               for (int v=0; v<NFLUX_TOTAL; v++)   CorrVal[v] = *FluidPtr1D[v] + FluxPtr[v][m][n]*Const[s];


//             calculate the internal energy density
#              if ( MODEL == HYDRO  &&  !defined BAROTROPIC_EOS && !defined SRHD )
               real Eint;
               real *ForEint = CorrVal;

//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
               real ForEint[NCOMP_TOTAL];
//             when FixSEint is on, use FluidPtr1D to calculate the original pressure
               if ( FixSEint )
                  for (int v=0; v<NCOMP_TOTAL; v++)   ForEint[v] = *FluidPtr1D[v];
               else
                  for (int v=0; v<NCOMP_TOTAL; v++)   ForEint[v] = CorrVal    [v];
*/

//             calculate the magnetic energy first
#              ifdef MHD
               int i, j, k;
               switch ( s )
               {
                  case 0:  i = 0;      j = n;      k = m;      break;
                  case 1:  i = PS1-1;  j = n;      k = m;      break;
                  case 2:  i = n;      j = 0;      k = m;      break;
                  case 3:  i = n;      j = PS1-1;  k = m;      break;
                  case 4:  i = n;      j = m;      k = 0;      break;
                  case 5:  i = n;      j = m;      k = PS1-1;  break;

                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
                     break;
               } // switch ( s )

               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
#              else
               const real Emag = NULL_REAL;
#              endif

//             when adopting the dual-energy formalism, we must determine to use Hydro_Con2Eint() or Hydro_DensEntropy2Pres()
//             since the fluid variables stored in CorrVal[] may not be fully consistent
//             --> because they have not been corrected by Hydro_DualEnergyFix()
//             --> also note that currently we adopt Hydro_DensEntropy2Pres() for DE_UPDATED_BY_MIN_PRES
//             --> consistency among all dual-energy related variables will be ensured after determining Eint
#              if ( DUAL_ENERGY == DE_ENPY )
               if ( *DE_StatusPtr1D == DE_UPDATED_BY_ETOT  ||  *DE_StatusPtr1D == DE_UPDATED_BY_ETOT_GRA )
#              endif
               {
                  const bool CheckMinEint_No = false;
                  Eint = Hydro_Con2Eint( ForEint[DENS], ForEint[MOMX], ForEint[MOMY], ForEint[MOMZ], ForEint[ENGY],
                                         NULL, NULL, CheckMinEint_No, NULL_REAL, Emag );
               }

#              if ( DUAL_ENERGY == DE_ENPY )
               else
               {
                  const bool CheckMinPres_No = false;
                  real Pres;
                  Pres = Hydro_DensEntropy2Pres( ForEint[DENS], ForEint[ENPY], EoS_AuxArray_Flt[1], CheckMinPres_No, NULL_REAL );
//                DE_ENPY only supports EOS_GAMMA, which does not involve passive scalars
                  Eint = EoS_DensPres2Eint_CPUPtr( ForEint[DENS], Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
               }
#              endif

#              if ( DUAL_ENERGY == DE_EINT )
#              error : DE_EINT is NOT supported yet !!
#              endif

#              endif // #if ( MODEL == HYDRO  &&  !defined BAROTROPIC_EOS )


//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
//             correct internal energy to restore the original specific internal energy
               if ( FixSEint )   Eint *= CorrVal[DENS] / *FluidPtr1D[DENS];
*/


//             do not apply the flux correction if there are any unphysical results
               bool ApplyFix;

#              if   ( MODEL == HYDRO )
#              ifdef SRHD
               if ( SRHD_CheckUnphysical( CorrVal, NULL,  __FUNCTION__, __LINE__, true  ) )
#              else
               if ( CorrVal[DENS] <= MIN_DENS
#                   ifndef BAROTROPIC_EOS
                    ||  Eint <= MIN_EINT  ||  !Aux_IsFinite(Eint)
#                   endif
#                   if   ( DUAL_ENERGY == DE_ENPY )
                    ||  ( (*DE_StatusPtr1D == DE_UPDATED_BY_DUAL || *DE_StatusPtr1D == DE_UPDATED_BY_MIN_PRES)
                           && CorrVal[ENPY] <= (real)2.0*TINY_NUMBER )

#                   elif ( DUAL_ENERGY == DE_EINT )
#                   error : DE_EINT is NOT supported yet !!
#                   endif
                  )
#              endif
#              elif ( MODEL == ELBDM  &&  defined CONSERVE_MASS )
               if ( CorrVal[DENS] <= MIN_DENS )
#              endif
               {
                  ApplyFix = false;
               }

               else
               {
                  ApplyFix = true;
               }


               if ( ApplyFix )
               {
//                floor and normalize the passive scalars
#                 if ( NCOMP_PASSIVE > 0 )
                  for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  CorrVal[v] = FMAX( CorrVal[v], TINY_NUMBER );

                  if ( OPT__NORMALIZE_PASSIVE )
                     Hydro_NormalizePassive( CorrVal[DENS], CorrVal+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#                 endif


//                ensure the consistency between pressure, total energy density, and dual-energy variable
//                --> assuming the variable "Eint" is correct
//                --> no need to check the internal energy floor here since we have skipped failing cells
#                 if ( MODEL == HYDRO && !defined SRHD )

//                for barotropic EoS, do not apply flux correction at all
#                 ifdef BAROTROPIC_EOS
                  CorrVal[ENGY] = *FluidPtr1D[ENGY];  // just set to the input value

#                 else
                  CorrVal[ENGY] = Hydro_ConEint2Etot( CorrVal[DENS], CorrVal[MOMX], CorrVal[MOMY], CorrVal[MOMZ], Eint, Emag );
#                 if   ( DUAL_ENERGY == DE_ENPY )
//                DE_ENPY only supports EOS_GAMMA, which does not involve passive scalars
                  CorrVal[ENPY] = Hydro_DensPres2Entropy( CorrVal[DENS],
                                                          EoS_DensEint2Pres_CPUPtr(CorrVal[DENS],Eint,NULL,
                                                          EoS_AuxArray_Flt,EoS_AuxArray_Int,h_EoS_Table),
                                                          EoS_AuxArray_Flt[1] );
#                 elif ( DUAL_ENERGY == DE_EINT )
#                 error : DE_EINT is NOT supported yet !!
#                 endif // DUAL_ENERGY

#                 endif // #ifdef BAROTROPIC_EOS ... else ...
#                 endif // #if ( MODEL == HYDRO )


//                store the corrected results
                  for (int v=0; v<NFLUX_TOTAL; v++)   *FluidPtr1D[v] = CorrVal[v];


//                rescale the real and imaginary parts to be consistent with the corrected amplitude
//                --> must NOT use CorrVal[REAL] and CorrVal[IMAG] below since NFLUX_TOTAL == 1 for ELBDM
#                 if ( MODEL == ELBDM  &&  defined CONSERVE_MASS )
                  real Re, Im, Rho_Corr, Rho_Wrong, Rescale;

                  Re        = *FluidPtr1D[REAL];
                  Im        = *FluidPtr1D[IMAG];
                  Rho_Corr  = *FluidPtr1D[DENS];
                  Rho_Wrong = SQR(Re) + SQR(Im);

//                be careful about the negative density introduced from the round-off errors
                  if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
                  {
                     *FluidPtr1D[DENS] = (real)0.0;
                     Rescale           = (real)0.0;
                  }
                  else
                     Rescale = SQRT( Rho_Corr/Rho_Wrong );

                  *FluidPtr1D[REAL] *= Rescale;
                  *FluidPtr1D[IMAG] *= Rescale;
#                 endif
               } // if ( ApplyFix )


//             update the fluid pointers
               for (int v=0; v<NCOMP_TOTAL; v++)
               FluidPtr1D[v]  += didx_n;

#              ifdef DUAL_ENERGY
               DE_StatusPtr1D += didx_n;
#              endif

            } // for (int n=0; n<PS1; n++}
         } // for (int m=0; m<PS1; m++}
      } // for (int s=0; s<6; s++)
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// 3. reset all flux arrays (in both real and buffer patches) to zero for bitwise reproducibility
#  ifdef BIT_REP_FLUX
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      for (int s=0; s<6; s++)
      {
         real (*FluxPtr)[PS1][PS1] = NULL;

         FluxPtr = amr->patch[0][lv][PID]->flux[s];
         if ( FluxPtr != NULL )
         {
            for (int v=0; v<NFLUX_TOTAL; v++)
            for (int m=0; m<PS1; m++)
            for (int n=0; n<PS1; n++)
               FluxPtr[v][m][n] = (real)0.0;
         }

         FluxPtr = amr->patch[0][lv][PID]->flux_bitrep[s];
         if ( FluxPtr != NULL )
         {
            for (int v=0; v<NFLUX_TOTAL; v++)
            for (int m=0; m<PS1; m++)
            for (int n=0; n<PS1; n++)
               FluxPtr[v][m][n] = (real)0.0;
         }
      }
   }
#  endif // #ifdef BIT_REP_FLUX

} // FUNCTION : Flu_FixUp_Flux
