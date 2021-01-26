#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Conservation
// Description :  Verify the conservation laws
//                --> HYDRO    : check mass, momentum, and energy
//                    ELBDM    : check mass and energy
//                    PAR_ONLY : check mass, momentum, and energy
//                    Passive scalars
//
// Note        :  1. This check only works with the models HYDRO, ELBDM, and PAR_ONLY
//                2. The values measured during the first function call will be taken as the reference values
//                   to estimate errors
//                   --> Note that during RESTART the reference values will be recalculated since they are NOT
//                       recorded in the output files currently
//                       --> Error estimation will be incorrect ...
//                3. For simulations with particles (i.e., when PARTICLE is on), the total conserved variables
//                   (e.g., total energy of gas and particles) will also be recorded
//
// Parameter   :  comment : You can put the location where this function is invoked in this string
//                          (not used currently)
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Conservation( const char *comment )
{

   static bool FirstTime = true;
   const char *FileName  = "Record__Conservation";


// check
#  ifdef COMOVING
   Aux_Message( stderr, "WARNING : function \"%s\" is NOT supported in COMOVING !!\n", __FUNCTION__ );
   OPT__CK_CONSERVATION = false;
   return;
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM  &&  MODEL != PAR_ONLY )
   Aux_Message( stderr, "WARNING : function \"%s\" is supported only in the models HYDRO, ELBDM, and PAR_ONLY !!\n",
                __FUNCTION__ );
   OPT__CK_CONSERVATION = false;
   return;
#  endif

   if ( FirstTime  &&  MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
   }


#  if   ( MODEL == HYDRO )
#  ifdef MHD
   const int    NVar_NoPassive    = 9;    // 9: mass, momentum (x/y/z), kinetic/internal/potential/magnetic/total energies
                                          // --> note that **total energy** is put in the last element
#  elif ( defined SRHD )
   const int    NVar_NoPassive    = 9;    // 8: mass, momentum (x/y/z), kinetic/internal/reduced/total energies

#  else                                   // --> note that **total energy** is put in the last element
   const int    NVar_NoPassive    = 8;    // 8: mass, momentum (x/y/z), kinetic/internal/potential/total energies
                                          // --> note that **total energy** is put in the last element
#  endif
   const int    idx_etot          = NVar_NoPassive - 1;
   const bool   CheckMinEint_No   = false;

#  elif ( MODEL == ELBDM )
   const int    NVar_NoPassive    = 5;    // 5: mass, kinetic/gravitational/self-interaction/total energies
   const int    idx_etot          = NVar_NoPassive - 1;
   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const int    NGhost            = 1;    // number of ghost zones for calculating the gradient of wave function
   const int    Size_Flu          = PS1 + 2*NGhost;
   const int    NPG               = 1;
   const double _2Eta2            = 0.5/SQR(ELBDM_ETA);
   const IntScheme_t IntScheme    = INT_CQUAR;

   real GradR[3], GradI[3], _dh2;
   int  ip, im, jp, jm, kp, km;

   real (*Flu_ELBDM)[2][Size_Flu][Size_Flu][Size_Flu] = new real [NPG*8][2][Size_Flu][Size_Flu][Size_Flu];

#  else
#  error : ERROR : unsupported MODEL !!
#  endif


// get the sum of passive scalars to be normalized
   const bool GetPassiveSum = ( PassiveNorm_NVar > 0 );
   const int  NVar_Max      = NVar_NoPassive + NCOMP_PASSIVE + 1; // for declaring the static variable Fluid_Ref
   const int  NVar          = NVar_NoPassive + NCOMP_PASSIVE + ( (GetPassiveSum)?1:0 );

   double dv, Fluid_ThisRank[NVar], Fluid_AllRank[NVar], Fluid_lv[NVar];   // dv : cell volume at each level
   int    FluSg;
#  ifdef GRAVITY
   int    PotSg;
#  endif
#  ifdef MHD
   int    MagSg;
#  endif
   FILE  *File = NULL;


// initialize accumulative variables as zero
   for (int v=0; v<NVar; v++)    Fluid_ThisRank[v] = 0.0;


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int v=0; v<NVar; v++)    Fluid_lv[v] = 0.0;

      dv    = CUBE( amr->dh[lv] );
      FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      PotSg = amr->PotSg[lv];
#     endif
#     ifdef MHD
      MagSg = amr->MagSg[lv];
#     endif
#     if ( MODEL == ELBDM )
      _dh2  = 0.5/amr->dh[lv];
#     endif


      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
#        if ( MODEL == ELBDM)
         const real MinDens_No = -1.0;
         const real MinPres_No = -1.0;

         Prepare_PatchData( lv, Time[lv], Flu_ELBDM[0][0][0][0], NULL, NGhost, NPG, &PID0, _REAL|_IMAG, _NONE,
                            IntScheme, INT_NONE, UNIT_PATCH, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );
#        endif

         for (int PID=PID0; PID<PID0+8; PID++)
         {
//          only check the leaf patches
            if ( amr->patch[0][lv][PID]->son == -1 )
            {
#              if   ( MODEL == HYDRO )
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  double Dens, MomX, MomY, MomZ, Etot, Ekin, Eint, Ered;
#                 ifdef GRAVITY
                  double Epot;
#                 endif
                  double Emag = NULL_REAL;

                  Dens = amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];
                  MomX = amr->patch[FluSg][lv][PID]->fluid[MOMX][k][j][i];
                  MomY = amr->patch[FluSg][lv][PID]->fluid[MOMY][k][j][i];
                  MomZ = amr->patch[FluSg][lv][PID]->fluid[MOMZ][k][j][i];
#                 ifdef SRHD
                  Ered = amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i];
#                 else
                  Etot = amr->patch[FluSg][lv][PID]->fluid[ENGY][k][j][i];
#                 endif

                  Fluid_lv[0] += Dens;
                  Fluid_lv[1] += MomX;
                  Fluid_lv[2] += MomY;
                  Fluid_lv[3] += MomZ;

#                 ifdef MHD
                  Emag         = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
                  Fluid_lv[7] += Emag;
#                 endif

#                 ifdef GRAVITY
//                set potential energy to zero when enabling both OPT__SELF_GRAVITY and OPT__EXT_POT
//                since the potential energy obtained here would be wrong anyway
//                --> to avoid possible misinterpretation
                  if      (  OPT__SELF_GRAVITY  &&  !OPT__EXT_POT )  Epot = 0.5*Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
                  else if ( !OPT__SELF_GRAVITY  &&   OPT__EXT_POT )  Epot =     Dens*amr->patch[PotSg][lv][PID]->pot[k][j][i];
                  else                                               Epot = 0.0;
                  Fluid_lv[6] += Epot;
#                 endif

#                 ifdef SRHD
                  printf("%s: %d need EoS table !! \n", __FUNCTION__, __LINE__);
                  Eint         = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Ered, EoS_GuessHTilde_CPUPtr,
                                                 EoS_HTilde2Temp_CPUPtr, CheckMinEint_No, NULL_REAL, Emag );
                  exit(0);
#                 else
                  Eint         = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Etot, EoS_GuessHTilde_CPUPtr,
                                                 EoS_HTilde2Temp_CPUPtr, CheckMinEint_No, NULL_REAL, Emag );
#                 endif

                  Fluid_lv[5] += Eint;

#                 ifdef SRHD
                  real Cons[NCOMP_FLUID] = { Dens, MomX, MomY, MomZ, Ered };
                  Ekin         = SRHD_Con2KineticEngy( Cons, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                  Fluid_lv[7] += Ered;
#                 else
                  //###NOTE: assuming Etot = Eint + Ekin + Emag
                  Ekin         = Etot - Eint;
#                 endif
#                 ifdef MHD
                  Ekin        -= Emag;
#                 endif
                  Fluid_lv[4] += Ekin;
               } // i,j,k


#              elif ( MODEL == ELBDM )
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
//                [0] mass
                  Fluid_lv[0] += amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i];

//                [2] potential energy in ELBDM
//                set potential energy to zero when enabling both OPT__SELF_GRAVITY and OPT__EXT_POT
//                since the potential energy obtained here would be wrong anyway
//                --> to avoid possible misinterpretation
#                 ifdef GRAVITY
                  if      (  OPT__SELF_GRAVITY  &&  !OPT__EXT_POT )
                  Fluid_lv[2] += 0.5*amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i]
                                    *amr->patch[PotSg][lv][PID]->pot        [k][j][i];

                  else if ( !OPT__SELF_GRAVITY  &&   OPT__EXT_POT )
                  Fluid_lv[2] +=     amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i]
                                    *amr->patch[PotSg][lv][PID]->pot        [k][j][i];
#                 endif

//                [3] quartic self-interaction potential
#                 ifdef QUARTIC_SELF_INTERACTION
                  Fluid_lv[3] += 0.5*ELBDM_LAMBDA*SQR( amr->patch[FluSg][lv][PID]->fluid[DENS][k][j][i] );
#                 endif
               }

//             [1] kinetic energy in ELBDM
               const int t = PID - PID0;

               for (int k=NGhost; k<Size_Flu-NGhost; k++)   { kp = k+1; km = k-1;
               for (int j=NGhost; j<Size_Flu-NGhost; j++)   { jp = j+1; jm = j-1;
               for (int i=NGhost; i<Size_Flu-NGhost; i++)   { ip = i+1; im = i-1;

                  GradR[0] = _dh2*( Flu_ELBDM[t][0][k ][j ][ip] - Flu_ELBDM[t][0][k ][j ][im] );
                  GradR[1] = _dh2*( Flu_ELBDM[t][0][k ][jp][i ] - Flu_ELBDM[t][0][k ][jm][i ] );
                  GradR[2] = _dh2*( Flu_ELBDM[t][0][kp][j ][i ] - Flu_ELBDM[t][0][km][j ][i ] );

                  GradI[0] = _dh2*( Flu_ELBDM[t][1][k ][j ][ip] - Flu_ELBDM[t][1][k ][j ][im] );
                  GradI[1] = _dh2*( Flu_ELBDM[t][1][k ][jp][i ] - Flu_ELBDM[t][1][k ][jm][i ] );
                  GradI[2] = _dh2*( Flu_ELBDM[t][1][kp][j ][i ] - Flu_ELBDM[t][1][km][j ][i ] );

                  Fluid_lv[1] += _2Eta2*( SQR(GradR[0]) + SQR(GradR[1]) + SQR(GradR[2]) +
                                          SQR(GradI[0]) + SQR(GradI[1]) + SQR(GradI[2])   );
               }}}


#              else
#              error : ERROR : unsupported MODEL !!

#              endif // MODEL


//             individual passive scalars
               for (int v=0; v<NCOMP_PASSIVE; v++)
               {
                  const int v1 = NVar_NoPassive + v;
                  const int v2 = NCOMP_FLUID    + v;

                  for (int k=0; k<PATCH_SIZE; k++)
                  for (int j=0; j<PATCH_SIZE; j++)
                  for (int i=0; i<PATCH_SIZE; i++)
                     Fluid_lv[v1] += amr->patch[FluSg][lv][PID]->fluid[v2][k][j][i];
               }
            } // if ( amr->patch[0][lv][PID]->son == -1 )
         } // for (int PID=PID0; PID<PID0+8; PID++)
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

//    get the total energy
#     if   ( MODEL == HYDRO )
#     ifdef MHD
      Fluid_lv[idx_etot] = Fluid_lv[4] + Fluid_lv[5] + Fluid_lv[6] + Fluid_lv[7];
#     elif ( defined SRHD )
      Fluid_lv[idx_etot] = Fluid_lv[4] + Fluid_lv[5] + Fluid_lv[6] + Fluid_lv[7];
#     else
      Fluid_lv[idx_etot] = Fluid_lv[4] + Fluid_lv[5] + Fluid_lv[6];
#     endif //#     ifdef MHD

#     elif ( MODEL == ELBDM )
      Fluid_lv[idx_etot] = Fluid_lv[1] + Fluid_lv[2] + Fluid_lv[3];
#     else
#     error : ERROR : unsupported MODEL !!
#     endif

//    sum of passive scalars to be normalized
      for (int v=0; v<PassiveNorm_NVar; v++)
         Fluid_lv[ NVar - 1 ] += Fluid_lv[ NVar_NoPassive + PassiveNorm_VarIdx[v] ];

//    sum over all levels
      for (int v=0; v<NVar; v++)    Fluid_ThisRank[v] += Fluid_lv[v]*dv;
   } // for (int lv=0; lv<NLEVEL; lv++)


// sum over all ranks
   MPI_Reduce( Fluid_ThisRank, Fluid_AllRank, NVar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// calculate conserved quantities for particles
#  ifdef PARTICLE
   double Mass_Par, MomX_Par, MomY_Par, MomZ_Par, Ekin_Par, Epot_Par, Etot_Par;

   Par_Aux_GetConservedQuantity( Mass_Par, MomX_Par, MomY_Par, MomZ_Par, Ekin_Par, Epot_Par );

   Etot_Par = Ekin_Par + Epot_Par;
#  endif


// output
   if ( MPI_Rank == 0 )
   {
//    note that a variable length array cannot have static storage duration
      static double Fluid_Ref[NVar_Max];
#     ifdef PARTICLE
      static double Mass_Par_Ref, MomX_Par_Ref, MomY_Par_Ref, MomZ_Par_Ref, Ekin_Par_Ref, Epot_Par_Ref, Etot_Par_Ref;
#     endif
      double AbsErr[NVar], RelErr[NVar];

      if ( FirstTime )
      {
//       record the reference values
         for (int v=0; v<NVar; v++)    Fluid_Ref[v] = Fluid_AllRank[v];

#        ifdef PARTICLE
         Mass_Par_Ref = Mass_Par;
         MomX_Par_Ref = MomX_Par;
         MomY_Par_Ref = MomY_Par;
         MomZ_Par_Ref = MomZ_Par;
         Ekin_Par_Ref = Ekin_Par;
         Epot_Par_Ref = Epot_Par;
         Etot_Par_Ref = Etot_Par;
#        endif


//       output header
         FILE *File = fopen( FileName, "a" );

         Aux_Message( File, "# Ref time     : %13.7e\n", Time[0] );
         Aux_Message( File, "# Ref step     : %ld\n",    Step    );
         Aux_Message( File, "\n" );

#        if   ( MODEL == HYDRO )
         Aux_Message( File, "# Mass_Gas     : total HYDRO mass\n" );
         Aux_Message( File, "# MomX/Y/Z_Gas : total HYDRO momentum\n" );
         Aux_Message( File, "# Ekin_Gas     : total HYDRO kinetic energy\n" );
         Aux_Message( File, "# Eint_Gas     : total HYDRO internal energy\n" );
         Aux_Message( File, "# Epot_Gas     : total HYDRO potential energy\n" );
#        ifdef MHD
         Aux_Message( File, "# Emag_Gas     : total HYDRO magnetic energy\n" );
#        endif
#        ifdef SRHD
         Aux_Message( File, "# Ered_Gas     : reduced HYDRO energy (rest mass energy subtracted from total energy)\n" );
#        endif
         Aux_Message( File, "# Etot_Gas     : total HYDRO energy\n" );

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "# Mass_Psi     : total ELBDM mass\n" );
         Aux_Message( File, "# Ekin_Psi     : total ELBDM kinetic energy\n" );
         Aux_Message( File, "# Epot_Psi     : total ELBDM potential energy\n" );
         Aux_Message( File, "# Esel_Psi     : total ELBDM self-interaction energy\n" );
         Aux_Message( File, "# Etot_Psi     : total ELBDM energy\n" );
#        endif

         if ( GetPassiveSum )
         Aux_Message( File, "# PassNorm     : sum of all target passive scalars to be normalized\n" );

#        ifdef PARTICLE
         Aux_Message( File, "# Mass_Par     : total PARTICLE mass\n" );
         Aux_Message( File, "# MomX/Y/Z_Par : total PARTICLE momentum\n" );
         Aux_Message( File, "# Ekin_Par     : total PARTICLE kinetic energy\n" );
         Aux_Message( File, "# Epot_Par     : total PARTICLE potential energy\n" );
         Aux_Message( File, "# Etot_Par     : total PARTICLE energy\n" );

#        if ( MODEL != PAR_ONLY )
         Aux_Message( File, "# Mass_All     : sum of the total HYDRO/ELBDM + PARTICLE mass\n" );
#        if ( MODEL == HYDRO )
         Aux_Message( File, "# MomX_All     : sum of the total HYDRO/ELBDM + PARTICLE momentum x\n" );
         Aux_Message( File, "# MomY_All     : sum of the total HYDRO/ELBDM + PARTICLE momentum y\n" );
         Aux_Message( File, "# MomZ_All     : sum of the total HYDRO/ELBDM + PARTICLE momentum z\n" );
#        endif
         Aux_Message( File, "# Etot_All     : sum of the total HYDRO/ELBDM + PARTICLE energy\n" );
#        endif // if ( MODEL != PAR_ONLY )
#        endif // #ifdef PARTICLE

         Aux_Message( File, "\n" );
         Aux_Message( File, "# AErr         : absolute error --> (now - ref)\n" );
         Aux_Message( File, "# RErr         : relative error --> (now - ref) / abs(ref)\n" );

         Aux_Message( File, "#-------------------------------------------------------------------------------------------" );
         Aux_Message( File, "--------------------------------------------------------------------------------------------\n\n" );

         const int NColumnMax = 80;
         Aux_Message( File, "#%12s  %10s", "[ 1]", "[ 2]" );
         for (int c=2; c<NColumnMax; c++)    Aux_Message( File, "  %10s[%2d]", "", c+1 );
         Aux_Message( File, "\n" );

         Aux_Message( File, "#%12s  %10s", "Time", "Step" );

#        if   ( MODEL == HYDRO )
         Aux_Message( File, "  %14s  %14s  %14s", "Mass_Gas", "Mass_Gas_AErr", "Mass_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomX_Gas", "MomX_Gas_AErr", "MomX_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomY_Gas", "MomY_Gas_AErr", "MomY_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomZ_Gas", "MomZ_Gas_AErr", "MomZ_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Ekin_Gas", "Ekin_Gas_AErr", "Ekin_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Eint_Gas", "Eint_Gas_AErr", "Eint_Gas_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Epot_Gas", "Epot_Gas_AErr", "Epot_Gas_RErr" );
#        ifdef MHD
         Aux_Message( File, "  %14s  %14s  %14s", "Emag_Gas", "Emag_Gas_AErr", "Emag_Gas_RErr" );
#        endif
#        ifdef SRHD
         Aux_Message( File, "  %14s  %14s  %14s", "Ered_Gas", "Ered_Gas_AErr", "Ered_Gas_RErr" );
#        endif
         Aux_Message( File, "  %14s  %14s  %14s", "Etot_Gas", "Etot_Gas_AErr", "Etot_Gas_RErr" );

#        elif ( MODEL == ELBDM )
         Aux_Message( File, "  %14s  %14s  %14s", "Mass_Psi", "Mass_Psi_AErr", "Mass_Psi_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Ekin_Psi", "Ekin_Psi_AErr", "Ekin_Psi_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Epot_Psi", "Epot_Psi_AErr", "Epot_Psi_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Esel_Psi", "Esel_Psi_AErr", "Esel_Psi_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Etot_Psi", "Etot_Psi_AErr", "Etot_Psi_RErr" );
#        endif

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         Aux_Message( File, "  %14s  %9s_AErr  %9s_RErr", FieldLabel[v], FieldLabel[v], FieldLabel[v] );

         if ( GetPassiveSum )
         Aux_Message( File, "  %14s  %14s  %14s", "PassNorm", "PassNorm_AErr", "PassNorm_RErr" );

#        ifdef PARTICLE
         Aux_Message( File, "  %14s  %14s  %14s", "Mass_Par", "Mass_Par_AErr", "Mass_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomX_Par", "MomX_Par_AErr", "MomX_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomY_Par", "MomY_Par_AErr", "MomY_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomZ_Par", "MomZ_Par_AErr", "MomZ_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Ekin_Par", "Ekin_Par_AErr", "Ekin_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Epot_Par", "Epot_Par_AErr", "Epot_Par_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "Etot_Par", "Etot_Par_AErr", "Etot_Par_RErr" );

#        if ( MODEL != PAR_ONLY )
         Aux_Message( File, "  %14s  %14s  %14s", "Mass_All", "Mass_All_AErr", "Mass_All_RErr" );
#        if ( MODEL == HYDRO )
         Aux_Message( File, "  %14s  %14s  %14s", "MomX_All", "MomX_All_AErr", "MomX_All_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomY_All", "MomY_All_AErr", "MomY_All_RErr" );
         Aux_Message( File, "  %14s  %14s  %14s", "MomZ_All", "MomZ_All_AErr", "MomZ_All_RErr" );
#        endif
         Aux_Message( File, "  %14s  %14s  %14s", "Etot_All", "Etot_All_AErr", "Etot_All_RErr" );
#        endif // if ( MODEL != PAR_ONLY )
#        endif // #ifdef PARTICLE

         Aux_Message( File, "\n" );

         fclose( File );
      } // if ( FirstTime )


//    calculate errors
      for (int v=0; v<NVar; v++)
      {
         AbsErr[v] = Fluid_AllRank[v] - Fluid_Ref[v];
         RelErr[v] = AbsErr[v] / fabs(Fluid_Ref[v]);
      }

//    calculate the sum of conserved quantities in different models
#     if ( defined PARTICLE  &&  MODEL != PAR_ONLY )
#     if   ( MODEL == HYDRO )
      const double Mass_All     = Fluid_AllRank[       0] + Mass_Par;
      const double MomX_All     = Fluid_AllRank[       1] + MomX_Par;
      const double MomY_All     = Fluid_AllRank[       2] + MomY_Par;
      const double MomZ_All     = Fluid_AllRank[       3] + MomZ_Par;
      const double Etot_All     = Fluid_AllRank[idx_etot] + Etot_Par;      // for HYDRO, total energy is stored in the last element

      const double Mass_All_Ref = Fluid_Ref    [       0] + Mass_Par_Ref;
      const double MomX_All_Ref = Fluid_Ref    [       1] + MomX_Par_Ref;
      const double MomY_All_Ref = Fluid_Ref    [       2] + MomY_Par_Ref;
      const double MomZ_All_Ref = Fluid_Ref    [       3] + MomZ_Par_Ref;
      const double Etot_All_Ref = Fluid_Ref    [idx_etot] + Etot_Par_Ref;

#     elif ( MODEL == ELBDM )
      const double Mass_All     = Fluid_AllRank[       0] + Mass_Par;
      const double Etot_All     = Fluid_AllRank[idx_etot] + Etot_Par;      // for ELBDM, total energy is stored in the last element

      const double Mass_All_Ref = Fluid_Ref    [       0] + Mass_Par_Ref;
      const double Etot_All_Ref = Fluid_Ref    [idx_etot] + Etot_Par_Ref;
#     endif // MODEL
#     endif // if ( defined PARTICLE  &&  MODEL != PAR_ONLY )


//    output
      File = fopen( FileName, "a" );

      Aux_Message( File, "%13.7e  %10ld", Time[0], Step );

      for (int v=0; v<NVar; v++)
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Fluid_AllRank[v], AbsErr[v], RelErr[v] );

#     ifdef PARTICLE
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Mass_Par, Mass_Par-Mass_Par_Ref, (Mass_Par-Mass_Par_Ref)/fabs(Mass_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomX_Par, MomX_Par-MomX_Par_Ref, (MomX_Par-MomX_Par_Ref)/fabs(MomX_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomY_Par, MomY_Par-MomY_Par_Ref, (MomY_Par-MomY_Par_Ref)/fabs(MomY_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomZ_Par, MomZ_Par-MomZ_Par_Ref, (MomZ_Par-MomZ_Par_Ref)/fabs(MomZ_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Ekin_Par, Ekin_Par-Ekin_Par_Ref, (Ekin_Par-Ekin_Par_Ref)/fabs(Ekin_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Epot_Par, Epot_Par-Epot_Par_Ref, (Epot_Par-Epot_Par_Ref)/fabs(Epot_Par_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Etot_Par, Etot_Par-Etot_Par_Ref, (Etot_Par-Etot_Par_Ref)/fabs(Etot_Par_Ref) );

#     if ( MODEL != PAR_ONLY )
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Mass_All, Mass_All-Mass_All_Ref, (Mass_All-Mass_All_Ref)/fabs(Mass_All_Ref) );
#     if ( MODEL == HYDRO )
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomX_All, MomX_All-MomX_All_Ref, (MomX_All-MomX_All_Ref)/fabs(MomX_All_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomY_All, MomY_All-MomY_All_Ref, (MomY_All-MomY_All_Ref)/fabs(MomY_All_Ref) );
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", MomZ_All, MomZ_All-MomZ_All_Ref, (MomZ_All-MomZ_All_Ref)/fabs(MomZ_All_Ref) );
#     endif
      Aux_Message( File, "  %14.7e  %14.7e  %14.7e", Etot_All, Etot_All-Etot_All_Ref, (Etot_All-Etot_All_Ref)/fabs(Etot_All_Ref) );
#     endif // if ( MODEL != PAR_ONLY )
#     endif // #ifdef PARTICLE

      Aux_Message( File, "\n" );

      fclose( File );
   } // if ( MPI_Rank == 0 )


   if ( FirstTime )  FirstTime = false;


#  if ( MODEL == ELBDM )
   if ( Flu_ELBDM != NULL )   delete [] Flu_ELBDM;
#  endif

} // FUNCTION : Aux_Check_Conservation
