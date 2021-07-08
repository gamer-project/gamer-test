#ifndef __CUFLU_FULLSTEPUPDATE__
#define __CUFLU_FULLSTEPUPDATE__



#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



// external functions
#ifdef __CUDACC__

#if ( NCOMP_PASSIVE > 0 )
# include "CUFLU_Shared_FluUtility.cu"
#endif

#ifdef DUAL_ENERGY
# include "CUFLU_Shared_DualEnergy.cu"
#endif

#endif // #ifdef __CUDACC__




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Note        :  1. This function is shared by MHM, MHM_RP, and CTU schemes
//                2. Invoke dual-energy check if DualEnergySwitch is on
//                3. If any unphysical fluid cells were found in patch, Hydro_FullStepUpdate will return instantly.
//
// Parameter   :  g_Input          : Array storing the input fluid data
//                g_Output         : Array to store the updated fluid data
//                g_DE_Status      : Array to store the dual-energy status
//                g_FC_B           : Array storing the updated face-centered B field
//                                   --> For the dual-energy formalism only
//                g_Flux           : Array storing the input face-centered fluxes
//                                   --> Accessed with the array stride N_FL_FLUX even thought its actually
//                                       allocated size is N_FC_FLUX^3
//                dt               : Time interval to advance solution
//                dh               : Cell size
//                MinDens/Eint     : Density and internal energy floors
//                DualEnergySwitch : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive      : true --> normalize passive scalars so that the sum of their mass density
//                                            is equal to the gas mass density
//                NNorm            : Number of passive scalars to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx          : Target variable indices to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_VarIdx"
//                EoS              : EoS object
//                                   --> Only for obtaining Gamma used by the dual-energy formalism
//                State            : (1/0) --> (Fail to update fluid patch/otherwise)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_FC_B[][ PS2P1*SQR(PS2) ], const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real dt, const real dh, const real MinDens, const real MinEint,
                           const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                           const EoS_t *EoS, char *State )
{

   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const real dt_dh        = dt/dh;

   real dFlux[3][NCOMP_TOTAL], Output_1Cell[NCOMP_TOTAL], Emag;


   const int size_ij = SQR(PS2);
   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
      const int i_out    = idx_out % PS2;
      const int j_out    = idx_out % size_ij / PS2;
      const int k_out    = idx_out / size_ij;

//    for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
#     ifdef MHD
      const int i_flux   = i_out + 1;
      const int j_flux   = j_out + 1;
      const int k_flux   = k_out + 1;
#     else
      const int i_flux   = i_out;
      const int j_flux   = j_out;
      const int k_flux   = k_out;
#     endif
      const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FL_FLUX, N_FL_FLUX );

      const int i_in     = i_out + FLU_GHOST_SIZE;
      const int j_in     = j_out + FLU_GHOST_SIZE;
      const int k_in     = k_out + FLU_GHOST_SIZE;
      const int idx_in   = IDX321( i_in, j_in, k_in, FLU_NXT, FLU_NXT );


//    1. calculate flux difference to update the fluid data
      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
#        ifdef MHD
         dFlux[d][v] = g_Flux[d][v][idx_flux] - g_Flux[d][v][ idx_flux - didx_flux[d] ];
#        else
         dFlux[d][v] = g_Flux[d][v][ idx_flux + didx_flux[d] ] - g_Flux[d][v][idx_flux];
#        endif
      }

      for (int v=0; v<NCOMP_TOTAL; v++)
         Output_1Cell[v] = g_Input[v][idx_in] - dt_dh*( dFlux[0][v] + dFlux[1][v] + dFlux[2][v] );


//    we no longer ensure positive density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//        because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
//    --> this consideration holds even when DUAL_ENERGY is adopted (e.g., when density is negative,
//        even when DUAL_ENERGY is on, we still want to try the 1st-order-flux correction before setting a floor value)
//    --> but for barotropic EoS, we apply Eint floor here to avoid any false alarm caused by Eint<0
#     ifdef BAROTROPIC_EOS
#     ifdef MHD
      Emag = MHD_GetCellCenteredBEnergy( g_FC_B[MAGX], g_FC_B[MAGY], g_FC_B[MAGZ],
                                         PS2, PS2, PS2, i_out, j_out, k_out );
#     else
      Emag = NULL_REAL;
#     endif
//    Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );
      Output_1Cell[ENGY] = Hydro_CheckMinEintInEngy( Output_1Cell[DENS], Output_1Cell[MOMX],
                                                     Output_1Cell[MOMY], Output_1Cell[MOMZ],
                                                     Output_1Cell[ENGY], MinEint, Emag );
#     endif // #ifdef BAROTROPIC_EOS


//    2. floor and normalize passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Output_1Cell[v] = FMAX( Output_1Cell[v], TINY_NUMBER );

      if ( NormPassive )
         Hydro_NormalizePassive( Output_1Cell[DENS], Output_1Cell+NCOMP_FLUID, NNorm, NormIdx );
#     endif


//    3. apply the dual-energy formalism to correct the internal energy
//    --> currently, even when UNSPLIT_GRAVITY is on (which would update the internal energy), we still invoke
//        Hydro_DualEnergyFix() here and will fix the internal energy in the gravity solver for cells updated
//        by the dual-energy formalism (i.e., for cells with their dual-energy status marked as DE_UPDATED_BY_DUAL)
//    --> this feature might be modified in the future
#     ifdef DUAL_ENERGY
//    B field must be updated in advance
#     ifdef MHD
      Emag = MHD_GetCellCenteredBEnergy( g_FC_B[MAGX], g_FC_B[MAGY], g_FC_B[MAGZ],
                                         PS2, PS2, PS2, i_out, j_out, k_out );
#     else
      Emag = NULL_REAL;
#     endif
//    we no longer apply density and pressure floors here since we want to enable 1st-order-flux correction for that
      const bool CheckMinPres_No = false;
//    Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );

      Hydro_DualEnergyFix( Output_1Cell[DENS], Output_1Cell[MOMX], Output_1Cell[MOMY], Output_1Cell[MOMZ],
                           Output_1Cell[ENGY], Output_1Cell[ENPY], g_DE_Status[idx_out],
                           EoS->AuxArrayDevPtr_Flt[1], EoS->AuxArrayDevPtr_Flt[2], CheckMinPres_No, NULL_REAL,
                           DualEnergySwitch, Emag );
#     endif // #ifdef DUAL_ENERGY


//    4. store results to the output array
      for (int v=0; v<NCOMP_TOTAL; v++)   g_Output[v][idx_out] = Output_1Cell[v];


//    5. check unphysical cells within a patch

      if ( State != NULL )
      {

        if(Hydro_CheckUnphysical( Output_1Cell, NULL, NULL, NULL, NULL, __FILE__, __FUNCTION__, __LINE__, false ))
        {
#         ifdef __CUDACC__
          atomicOr ( (int*)State, 1);
#         else
          *State = *State | 1;
#         endif
        }


//      5-1. waiting all threads within a GPU block
#       ifdef __CUDACC__
        __syncthreads();
#       endif


//      5-2. return all threads within a block
        if ( *State == 1 )     return;

      }


//    6. check the negative density and energy
#     ifdef CHECK_UNPHYSICAL_IN_FLUID
      Hydro_CheckUnphysical( NULL, NULL, &Output_1Cell[DENS], NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
      Hydro_CheckUnphysical( NULL, NULL, &Output_1Cell[ENGY], NULL, "energy", __FILE__, __FUNCTION__, __LINE__, true );
#     endif

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

} // FUNCTION : Hydro_FullStepUpdate



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )



#endif // #ifndef __CUFLU_FULLSTEPUPDATE__
