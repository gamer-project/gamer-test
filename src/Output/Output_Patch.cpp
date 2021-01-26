#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Patch
// Description :  Output the data of a single patch
//
// Example     :  const int lv  = 2;
//                const int PID = 20;
//                char comment[10];
//                sprintf( comment, "Step%ld", AdvanceCounter[lv] );
//                Output_Patch( lv, PID, amr->FluSg[lv], amr->MagSg[lv], amr->PotSg[lv], comment );
//
// Parameter   :  lv      : Target refinement level
//                PID     : Target patch index
//                FluSg   : Sandglass of the fluid data
//                MagSg   : Sandglass of the magnetic field data
//                PotSg   : Sandglass of the potential data
//                comment : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Patch( const int lv, const int PID, const int FluSg, const int MagSg, const int PotSg, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( PID < 0  ||  PID >= MAX_PATCH )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (MAX_PATCH = %d) !!\n", "PID", PID, MAX_PATCH );

   if ( FluSg < 0  ||  FluSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  ifdef MHD
   if ( MagSg < 0  ||  MagSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "MagSg", MagSg );
#  endif

#  ifdef GRAVITY
   if ( PotSg < 0  ||  PotSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );
#  endif

   if ( amr->patch[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : level = %d, PID = %d does NOT exist !!\n", lv, PID );
      return;
   }


   patch_t *Relation                  = amr->patch[    0][lv][PID];
   real    (*fluid)[PS1][PS1][PS1]    = amr->patch[FluSg][lv][PID]->fluid;
#  ifdef MHD
   real    (*magnetic)[PS1P1*PS1*PS1] = amr->patch[MagSg][lv][PID]->magnetic;
#  endif
#  ifdef GRAVITY
   real    (*pot)[PS1][PS1]           = amr->patch[PotSg][lv][PID]->pot;
#  endif

   char FileName[100];
   sprintf( FileName, "Patch_r%d_lv%d_p%d", MPI_Rank, lv, PID );
   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// 1. output patch information
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  FluSg %d  PotSg %d  Time %13.7e  Step %ld  Counter %ld\n",
            MPI_Rank, lv, PID, PID%8, FluSg, PotSg, Time[lv], Step, AdvanceCounter[lv] );

   fprintf( File, "Father %d  Son %d  Corner (%10d,%10d,%10d)  Size %13.7e", Relation->father, Relation->son,
            Relation->corner[0], Relation->corner[1], Relation->corner[2], PS1*amr->dh[lv] );
#  ifdef LOAD_BALANCE
   fprintf( File, "  LB_Idx %ld  PaddedCr1D %lu", Relation->LB_Idx, Relation->PaddedCr1D );
#  endif
   fprintf( File, "\n" );
   fprintf( File, "EdgeL = (%20.14e, %20.14e, %20.14e)\n", Relation->EdgeL[0], Relation->EdgeL[1], Relation->EdgeL[2] );
   fprintf( File, "EdgeR = (%20.14e, %20.14e, %20.14e)\n", Relation->EdgeR[0], Relation->EdgeR[1], Relation->EdgeR[2] );

#  ifdef PARTICLE
   fprintf( File, "NPar = %5d  ParListSize = %5d\n", Relation->NPar, Relation->ParListSize );
#  endif


   fprintf( File, "\nSibling, Sibling->Son, and Father->Sibling Lists :\n" );

   int Sib, FaSib, SibSon, Fa;
   for (int S=0; S<26; S++)
   {
      Fa     = Relation->father;
      Sib    = Relation->sibling[S];
      FaSib  = ( Fa == -1 ) ? -1 : ( amr->patch[0][lv-1][Fa] != NULL ) ?
                                     amr->patch[0][lv-1][Fa]->sibling[S] : -999;
      SibSon = ( Sib < 0 )  ? Sib : amr->patch[0][lv][Sib]->son;

      fprintf( File, "Sib[%2d] = %6d     Sib_Son = %6d     Fa_Sib[%2d] = %6d\n",
               S, Sib, SibSon, S, FaSib );
   }
   fprintf( File, "\n" );



// check whether or not the target patch stores physical data
   if ( fluid == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store fluid data !!\n", lv, PID );
#  ifdef MHD
   if ( magnetic == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store magnetic field data !!\n", lv, PID );
#  endif
#  ifdef GRAVITY
   if ( pot == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store potential data !!\n", lv, PID );
#  endif



// 2. output physical data
// header
   fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );

   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, "%14s", FieldLabel[v] );

#  if   ( MODEL == HYDRO )
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)     fprintf( File, "%14s", MagLabel[v] );
   fprintf( File, "%14s", "MagEngy" );
#  endif
   fprintf( File, "%14s", "Pressure" );
#  ifdef DUAL_ENERGY
   fprintf( File, "%14s", "DE-status" );
#  endif

#  elif ( MODEL == ELBDM )

#  else
#  warning : WARNING : DO YOU WANT TO ADD THE FILE HEADER HERE FOR THE NEW MODEL ??
#  endif // MODEL

#  ifdef GRAVITY
   fprintf( File, "%14s", PotLabel );
#  endif

   fprintf( File, "\n" );


// output data
   real u[NCOMP_TOTAL];

   for (int k=0; k<PS1; k++)
   for (int j=0; j<PS1; j++)
   for (int i=0; i<PS1; i++)
   {
//    cell indices
      fprintf( File, "(%2d,%2d,%2d)", i, j, k );

      if ( fluid != NULL )
      {
//       output all variables in the fluid array
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            u[v] = fluid[v][k][j][i];
            fprintf( File, " %13.6e", u[v] );
         }

#        if   ( MODEL == HYDRO )
//       magnetic field
#        ifdef MHD
         const real Emag = ( magnetic == NULL ) ?
                           NULL_REAL :
                           MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg );
         real B[3] = { NULL_REAL, NULL_REAL, NULL_REAL };
         if ( magnetic != NULL )    MHD_GetCellCenteredBFieldInPatch( B, lv, PID, i, j, k, MagSg );
         fprintf( File, " %13.6e %13.6e %13.6e %13.6e", B[MAGX], B[MAGY], B[MAGZ], Emag );
#        endif

//       pressure
         const bool CheckMinPres_No = false;
#        ifdef MHD
//       set pressure to NULL_REAL if somehow magnetic[] is not allocated (likely due to a bug)
         const real Pres = ( magnetic == NULL ) ?
                           NULL_REAL :
                           Hydro_Con2Pres( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                                           CheckMinPres_No, NULL_REAL, Emag,
                                           EoS_DensEint2Pres_CPUPtr,
                                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
#        else
         const real Pres = Hydro_Con2Pres( u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], u+NCOMP_FLUID,
                                           CheckMinPres_No, NULL_REAL, NULL_REAL,
                                           EoS_DensEint2Pres_CPUPtr,
                                           EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                           EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
#        endif
         fprintf( File, " %13.6e", Pres );

//       dual-energy variable
#        ifdef DUAL_ENERGY
         fprintf( File, " %13c", Relation->de_status[k][j][i] );
#        endif

#        elif ( MODEL == ELBDM )

#        else
#        warning : WARNING : DO YOU WANT TO ADD THE MODEL-SPECIFIC FIELD HERE FOR THE NEW MODEL ??
#        endif // MODEL
      } // if ( fluid != NULL )

      else
      {
//       output empty strings if the fluid array is not allocated
         for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13s", "" );

#        if ( MODEL == HYDRO )
#        ifdef MHD
         fprintf( File, " %13s %13s %13s %13s", "", "", "", "" );
#        endif
         fprintf( File, " %13s", "" );
#        ifdef DUAL_ENERGY
         fprintf( File, " %13s", "" );
#        endif

#        elif ( MODEL == ELBDM )

#        else
#        warning : WARNING : DO YOU WANT TO ADD THE MODEL-SPECIFIC FIELD HERE FOR THE NEW MODEL ??
#        endif // MODEL
      } // if ( fluid != NULL ) ... else ...

//    potential
#     ifdef GRAVITY
      if ( pot != NULL )   fprintf( File, " %13.6e", pot[k][j][i] );
      else                 fprintf( File, " %13s", "" );
#     endif

      fprintf( File, "\n" );
   } // i,j,k



// 3. output face-centered magnetic field
#  ifdef MHD
   fprintf( File, "\n" );
   fprintf( File, "====================================\n" );
   fprintf( File, "== MAGNETIC FIELD (face-centered) == \n" );
   fprintf( File, "====================================\n" );
   fprintf( File, "\n" );

   if ( magnetic != NULL )
   {
//    header
      fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s\n", "i", "j", "k", MagLabel[MAGX], MagLabel[MAGY], MagLabel[MAGZ] );

      for (int k=0; k<PS1P1; k++)
      for (int j=0; j<PS1P1; j++)
      for (int i=0; i<PS1P1; i++)
      {
//       cell indices
         fprintf( File, "(%2d,%2d,%2d)", i, j, k );

//       B_X
         if ( j < PS1  &&  k < PS1 )   fprintf( File, " %13.6e", magnetic[MAGX][ IDX321_BX(i,j,k,PS1,PS1) ] );
         else                          fprintf( File, " %13s", "" );

//       B_Y
         if ( i < PS1  &&  k < PS1 )   fprintf( File, " %13.6e", magnetic[MAGY][ IDX321_BY(i,j,k,PS1,PS1) ] );
         else                          fprintf( File, " %13s", "" );

//       B_Z
         if ( i < PS1  &&  j < PS1 )   fprintf( File, " %13.6e", magnetic[MAGZ][ IDX321_BZ(i,j,k,PS1,PS1) ] );
         else                          fprintf( File, " %13s", "" );

         fprintf( File, "\n" );
      } // i,j,k
   } // if ( magnetic != NULL )
#  endif // #ifdef MHD



// 4. output particles
#  ifdef PARTICLE
   long ParID;
   fprintf( File, "\n" );
   fprintf( File, "===================\n" );
   fprintf( File, "== PARTICLE DATA == \n" );
   fprintf( File, "===================\n" );
   fprintf( File, "\n" );
   fprintf( File, "%5s  %10s", "No.", "ParID" );
   for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %13s", ParAttLabel[v] );
   fprintf( File, "\n" );

   for (int p=0; p<Relation->NPar; p++)
   {
      ParID = Relation->ParList[p];

      fprintf( File, "%5d  %10ld", p, ParID );
      for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %13.6e", amr->Par->Attribute[v][ParID] );

      fprintf( File, "\n" );
   }
#  endif // #ifdef PARTICLE


   fclose( File );

} // FUNCTION : Output_Patch

