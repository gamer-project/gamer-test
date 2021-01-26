#include "GAMER.h"

extern void SetTempIntPara( const int lv, const int Sg_Current, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeProfile
// Description :  Compute the average radial profile of target field(s)
//
// Note        :  1. Results will be stored in the input "Prof" object
//                   --> Prof->Radius[]: Radial coordinate at each bin
//                       Prof->Data  []: Profile data at each bin
//                       Prof->Weight[]: Total weighting at each bin
//                       Prof->NCell []: Number of cells at each bin
//                       Prof->NBin    : Total number of bins
//                   --> See the "Profile_t" structure defined in "include/Profile.h" for details
//                   --> These arrays will be free'd when deleting "Prof"
//                2. Maximum radius adopted when actually computing the profile may be larger than the input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//                3. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same profile data after invoking this function
//                4. Use cell volume as the weighting of each cell
//                   --> Will support other weighting functions in the future
//                5. Support computing multiple fields
//                   --> The order of fields to be returned follows TVarBitIdx[]
//
// Parameter   :  Prof        : Profile_t object array to store the results
//                Center      : Target center coordinates
//                r_max_input : Maximum radius for computing the profile
//                              --> See also "Note-2" above
//                dr_min      : Minimum bin size
//                              --> For linear bin, this is the size of all bins
//                                  For log    bin, this is the size of the 0th bin
//                LogBin      : true/false --> log/linear bins
//                LogBinRatio : Ratio of adjacent log bins
//                              --> Right edge of log bin n = dr_min*LogBinRatio^n
//                RemoveEmpty : true  --> remove empty bins from the data
//                              false --> these empty bins will still be in the profile arrays with
//                                        Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0
//                TVarBitIdx  : Bitwise indices of target variables for computing the profiles
//                              --> Supported indices (defined in Macro.h):
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELR, _PRES, _EINT_DER
//                                             [, _ENPY, _EINT, _POTE]
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                              --> For a passive scalar with an integer field index FieldIdx returned by AddField(),
//                                  one can convert it to a bitwise field index by BIDX(FieldIdx)
//                NProf       : Number of Profile_t objects in Prof
//                SingleLv    : Only consider patches on the specified level
//                              --> If SingleLv<0, loop over all levels
//                MaxLv       : Consider patches on levels equal/below MaxLv if SingleLv<0
//                              --> If MaxLv<0, loop over all levels
//                PatchType   : Only consider patches of the specified type
//                              --> Supported types: PATCH_LEAF, PATCH_NONLEAF, PATCH_BOTH
//                PrepTime    : Target physical time to prepare data
//                              --> If PrepTime<0, turn off temporal interpolation and always use the most recent data
//
// Example     :  const double      Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double      MaxRadius      = 0.5*amr->BoxSize[0];
//                const double      MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool        LogBin         = true;
//                const double      LogBinRatio    = 1.25;
//                const bool        RemoveEmptyBin = true;
//                const long        TVar[]         = { _DENS, _PRES };
//                const int         NProf          = 2;
//                const int         SingleLv       = -1;
//                const int         MaxLv          = -1;
//                const PatchType_t PatchType      = PATCH_LEAF;
//                const double      PrepTime       = -1.0;
//
//                Profile_t Prof_Dens, Prof_Pres;
//                Profile_t *Prof[] = { &Prof_Dens, &Prof_Pres };
//
//                Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
//                                    TVar, NProf, SingleLv, MaxLv, PatchType, PrepTime );
//
//                if ( MPI_Rank == 0 )
//                {
//                   for (int p=0; p<NProf; p++)
//                   {
//                      char Filename[MAX_STRING];
//                      sprintf( Filename, "Profile%d.txt", p );
//                      FILE *File = fopen( Filename, "w" );
//                      fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
//                      for (int b=0; b<Prof[p]->NBin; b++)
//                         fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
//                                  Prof[p]->Radius[b], Prof[p]->Data[b], Prof[p]->Weight[b], Prof[p]->NCell[b] );
//                      fclose( File );
//                   }
//                }
//
// Return      :  Prof
//-------------------------------------------------------------------------------------------------------
void Aux_ComputeProfile( Profile_t *Prof[], const double Center[], const double r_max_input, const double dr_min,
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const long TVarBitIdx[],
                         const int NProf, const int SingleLv, const int MaxLv, const PatchType_t PatchType,
                         const double PrepTime )
{

// check
#  ifdef GAMER_DEBUG
   if ( r_max_input <= 0.0 )
      Aux_Error( ERROR_INFO, "r_max_input (%14.7e) <= 0.0 !!\n", r_max_input );

   if ( dr_min <= 0.0 )
      Aux_Error( ERROR_INFO, "dr_min (%14.7e) <= 0.0 !!\n", dr_min );

   if ( LogBin  &&  LogBinRatio <= 1.0 )
      Aux_Error( ERROR_INFO, "LogBinRatio (%14.7e) <= 1.0 !!\n", LogBinRatio );

   if ( ( SingleLv >= 0 )  &&  ( MaxLv >= 0 ) )
      Aux_Error( ERROR_INFO, "SingleLv (%d) and MaxLv (%d) cannot be both >= 0 !!\n", SingleLv, MaxLv );
#  endif


// precompute the integer indices of intrinsic fluid fields for better performance
   const int IdxUndef = -1;
   int TFluIntIdx[NProf];

   for (int p=0; p<NProf; p++)
   {
      TFluIntIdx[p] = IdxUndef;

      for (int v=0; v<NCOMP_TOTAL; v++)
         if ( TVarBitIdx[p] & (1L<<v) )   TFluIntIdx[p] = v;
   }


// check whether _POTE is in TVarBitIdx since the potential array may have not been computed during initialization
#  ifdef GRAVITY
   bool InclPot = false;

   for (int p=0; p<NProf; p++)
      if ( TVarBitIdx[p] & _POTE )   InclPot = true;
#  endif


// initialize the profile objects
   for (int p=0; p<NProf; p++)
   {
//    get the total number of radial bins and the corresponding maximum radius
      if ( LogBin )
      {
         Prof[p]->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
         Prof[p]->MaxRadius = dr_min*pow( LogBinRatio, Prof[p]->NBin-1 );
      }

      else // linear bin
      {
         Prof[p]->NBin      = (int)ceil( r_max_input / dr_min );
         Prof[p]->MaxRadius = dr_min*Prof[p]->NBin;
      }


//    record profile parameters
      for (int d=0; d<3; d++)    Prof[p]->Center[d] = Center[d];

      Prof[p]->LogBin = LogBin;

      if ( LogBin )  Prof[p]->LogBinRatio = LogBinRatio;


//    allocate all member arrays of Prof
      Prof[p]->AllocateMemory();


//    record radial coordinates
      if ( LogBin )
         for (int b=0; b<Prof[0]->NBin; b++)    Prof[p]->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
      else
         for (int b=0; b<Prof[0]->NBin; b++)    Prof[p]->Radius[b] = (b+0.5)*dr_min;

   } // for (int p=0; p<NProf; p++)


// allocate memory for the per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double ***OMP_Data=NULL, ***OMP_Weight=NULL;
   long   ***OMP_NCell=NULL;

   Aux_AllocateArray3D( OMP_Data,   NProf, NT, Prof[0]->NBin );
   Aux_AllocateArray3D( OMP_Weight, NProf, NT, Prof[0]->NBin );
   Aux_AllocateArray3D( OMP_NCell,  NProf, NT, Prof[0]->NBin );


// collect profile data in this rank
   const double r_max2      = SQR( Prof[0]->MaxRadius );
   const double HalfBox[3]  = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[4] == BC_FLU_PERIODIC };

#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         OMP_Data  [p][TID][b] = 0.0;
         OMP_Weight[p][TID][b] = 0.0;
         OMP_NCell [p][TID][b] = 0;
      }

//    allocate passive scalar arrays
#     if ( MODEL == HYDRO )
      real *Passive      = new real [NCOMP_PASSIVE];
      real *Passive_IntT = new real [NCOMP_PASSIVE];
#     endif

//    determine which levels to be considered
      const int lv_min = ( SingleLv < 0 ) ? 0                                     : SingleLv;
      const int lv_max = ( SingleLv < 0 ) ? ( ( MaxLv < 0 ) ? TOP_LEVEL : MaxLv ) : SingleLv;

      for (int lv=lv_min; lv<=lv_max; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );


//       determine temporal interpolation parameters
         bool FluIntTime = false;
         int  FluSg      = amr->FluSg[lv];
         int  FluSg_IntT;
         real FluWeighting, FluWeighting_IntT;

#        ifdef MHD
         bool MagIntTime = false;
         int  MagSg      = amr->MagSg[lv];
         int  MagSg_IntT;
         real MagWeighting, MagWeighting_IntT;
#        endif

#        ifdef GRAVITY
         bool PotIntTime = false;
         int  PotSg      = amr->PotSg[lv];
         int  PotSg_IntT;
         real PotWeighting, PotWeighting_IntT;
#        endif

         if ( PrepTime >= 0.0 )
         {
//          fluid
            SetTempIntPara( lv, amr->FluSg[lv], PrepTime, amr->FluSgTime[lv][0], amr->FluSgTime[lv][1],
                            FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );

//          magnetic field
#           ifdef MHD
            SetTempIntPara( lv, amr->MagSg[lv], PrepTime, amr->MagSgTime[lv][0], amr->MagSgTime[lv][1],
                            MagIntTime, MagSg, MagSg_IntT, MagWeighting, MagWeighting_IntT );
#           endif

//          potential
#           ifdef GRAVITY
            if ( InclPot )
               SetTempIntPara( lv, amr->PotSg[lv], PrepTime, amr->PotSgTime[lv][0], amr->PotSgTime[lv][1],
                               PotIntTime, PotSg, PotSg_IntT, PotWeighting, PotWeighting_IntT );
#           endif
         }


//       use the "static" schedule for reproducibility
#        pragma omp for schedule( static )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
//          determine which type of patches to be looped
            if (  ( amr->patch[0][lv][PID]->son != -1 && PatchType == PATCH_LEAF    )  ||
                  ( amr->patch[0][lv][PID]->son == -1 && PatchType == PATCH_NONLEAF )  )
               continue;


            const real (*FluidPtr)[PS1][PS1][PS1] = amr->patch[ FluSg ][lv][PID]->fluid;
#           ifdef GRAVITY
            const real (*PotPtr  )[PS1][PS1]      = amr->patch[ PotSg ][lv][PID]->pot;
#           endif

//          pointer for temporal interpolation
            const real (*FluidPtr_IntT)[PS1][PS1][PS1] = ( FluIntTime ) ? amr->patch[ FluSg_IntT ][lv][PID]->fluid : NULL;
#           ifdef GRAVITY
            const real (*PotPtr_IntT  )[PS1][PS1]      = ( PotIntTime ) ? amr->patch[ PotSg_IntT ][lv][PID]->pot   : NULL;
#           endif


            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

            for (int k=0; k<PS1; k++)  {  double dz = z0 + k*dh;
                                          if ( Periodic[2] ) {
                                             if      ( dz > +HalfBox[2] )  {  dz -= amr->BoxSize[2];  }
                                             else if ( dz < -HalfBox[2] )  {  dz += amr->BoxSize[2];  }
                                          }
            for (int j=0; j<PS1; j++)  {  double dy = y0 + j*dh;
                                          if ( Periodic[1] ) {
                                             if      ( dy > +HalfBox[1] )  {  dy -= amr->BoxSize[1];  }
                                             else if ( dy < -HalfBox[1] )  {  dy += amr->BoxSize[1];  }
                                          }
            for (int i=0; i<PS1; i++)  {  double dx = x0 + i*dh;
                                          if ( Periodic[0] ) {
                                             if      ( dx > +HalfBox[0] )  {  dx -= amr->BoxSize[0];  }
                                             else if ( dx < -HalfBox[0] )  {  dx += amr->BoxSize[0];  }
                                          }

               const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

               if ( r2 < r_max2 )
               {
                  const double r   = sqrt( r2 );
                  const int    bin = ( LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(LogBinRatio) ) + 1  )
                                                : int( r/dr_min );
//                prevent from round-off errors
                  if ( bin >= Prof[0]->NBin )   continue;

//                check
#                 ifdef GAMER_DEBUG
                  if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                 endif

//                prepare passive scalars (for better sustainability, always do it even when unnecessary)
#                 if ( MODEL == HYDRO )
                  for (int v_out=0; v_out<NCOMP_PASSIVE; v_out++)
                  {
                     const int v_in = v_out + NCOMP_FLUID;

                     Passive     [v_out] = FluidPtr     [v_in][k][j][i];
                     if ( FluIntTime )
                     Passive_IntT[v_out] = FluidPtr_IntT[v_in][k][j][i];
                  }
#                 endif

                  for (int p=0; p<NProf; p++)
                  {
//                   intrinsic fluid fields
                     if ( TFluIntIdx[p] != IdxUndef )
                     {
                        const real Weight = dv;

                        OMP_Data  [p][TID][bin] += ( FluIntTime )
                                                 ? ( FluWeighting     *FluidPtr     [ TFluIntIdx[p] ][k][j][i]
                                                   + FluWeighting_IntT*FluidPtr_IntT[ TFluIntIdx[p] ][k][j][i] )*Weight
                                                 :                     FluidPtr     [ TFluIntIdx[p] ][k][j][i]  *Weight;
                        OMP_Weight[p][TID][bin] += Weight;
                        OMP_NCell [p][TID][bin] ++;
                     }

//                   other fields
                     else
                     {
                        switch ( TVarBitIdx[p] )
                        {
//                         gravitational potential
#                          ifdef GRAVITY
                           case _POTE:
                           {
                              const real Weight = ( FluIntTime )    // weighted by cell mass
                                                ? ( FluWeighting     *FluidPtr     [DENS][k][j][i]
                                                  + FluWeighting_IntT*FluidPtr_IntT[DENS][k][j][i] )*dv
                                                :                     FluidPtr     [DENS][k][j][i]  *dv;

                              OMP_Data  [p][TID][bin] += ( PotIntTime )
                                                       ? ( PotWeighting     *PotPtr     [k][j][i]
                                                         + PotWeighting_IntT*PotPtr_IntT[k][j][i] )*Weight
                                                       :                     PotPtr     [k][j][i]  *Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;
#                          endif

//                         derived fields
#                          if ( MODEL == HYDRO )
                           case _VELR:
                           {
                              const real Weight = ( FluIntTime )    // weighted by cell mass
                                                ? ( FluWeighting     *FluidPtr     [DENS][k][j][i]
                                                  + FluWeighting_IntT*FluidPtr_IntT[DENS][k][j][i] )*dv
                                                :                     FluidPtr     [DENS][k][j][i]  *dv;

                              const real MomR   = ( FluIntTime )
                                                ? ( FluWeighting     *( FluidPtr     [MOMX][k][j][i]*dx +
                                                                        FluidPtr     [MOMY][k][j][i]*dy +
                                                                        FluidPtr     [MOMZ][k][j][i]*dz )
                                                  + FluWeighting_IntT*( FluidPtr_IntT[MOMX][k][j][i]*dx +
                                                                        FluidPtr_IntT[MOMY][k][j][i]*dy +
                                                                        FluidPtr_IntT[MOMZ][k][j][i]*dz ) ) / r
                                                :                     ( FluidPtr     [MOMX][k][j][i]*dx +
                                                                        FluidPtr     [MOMY][k][j][i]*dy +
                                                                        FluidPtr     [MOMZ][k][j][i]*dz )   / r;

                              OMP_Data  [p][TID][bin] += MomR*dv;    // vr*(rho*dv)
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;

                           case _PRES:
                           {
                              const bool CheckMinPres_No = false;
                              const real Weight          = dv;
#                             ifdef MHD
                              const real Emag            = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg      );
                              const real Emag_IntT       = ( MagIntTime )
                                                         ? MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT )
                                                         : NULL_REAL;
#                             else
                              const real Emag            = NULL_REAL;
                              const real Emag_IntT       = NULL_REAL;
#                             endif
                              const real Pres = ( FluIntTime )
                                              ?   FluWeighting     *Hydro_Con2Pres( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    Passive,
                                                                                    CheckMinPres_No, NULL_REAL, Emag,
                                                                                    EoS_DensEint2Pres_CPUPtr,
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL )
                                                + FluWeighting_IntT*Hydro_Con2Pres( FluidPtr_IntT[DENS][k][j][i],
                                                                                    FluidPtr_IntT[MOMX][k][j][i],
                                                                                    FluidPtr_IntT[MOMY][k][j][i],
                                                                                    FluidPtr_IntT[MOMZ][k][j][i],
                                                                                    FluidPtr_IntT[ENGY][k][j][i],
                                                                                    Passive_IntT,
                                                                                    CheckMinPres_No, NULL_REAL, Emag_IntT,
                                                                                    EoS_DensEint2Pres_CPUPtr,
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL )
                                              :                     Hydro_Con2Pres( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    Passive,
                                                                                    CheckMinPres_No, NULL_REAL, Emag,
                                                                                    EoS_DensEint2Pres_CPUPtr,
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    EoS_AuxArray_Flt,
                                                                                    EoS_AuxArray_Int, h_EoS_Table, NULL );

                              OMP_Data  [p][TID][bin] += Pres*Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;

                           case _EINT_DER:
                           {
                              const real Weight = dv;
                              const real Dens   = FluidPtr[DENS][k][j][i];

//                            use the dual-energy variable to calculate the internal energy directly, if applicable
#                             ifdef DUAL_ENERGY

#                             if   ( DUAL_ENERGY == DE_ENPY )
                              const bool CheckMinPres_No = false;
                              const real Enpy = FluidPtr[ENPY][k][j][i];
                              const real Pres = Hydro_DensEntropy2Pres( Dens, Enpy, EoS_AuxArray_Flt[1],
                                                                        CheckMinPres_No, NULL_REAL );
                              const real Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, Passive, EoS_AuxArray_Flt,
                                                                          EoS_AuxArray_Int, h_EoS_Table );
#                             elif ( DUAL_ENERGY == DE_EINT )
#                             error : DE_EINT is NOT supported yet !!
#                             endif

#                             else // #ifdef DUAL_ENERGY

                              const bool CheckMinEint_No = false;
                              const real MomX            = FluidPtr[MOMX][k][j][i];
                              const real MomY            = FluidPtr[MOMY][k][j][i];
                              const real MomZ            = FluidPtr[MOMZ][k][j][i];
                              const real Etot            = FluidPtr[ENGY][k][j][i];
#                             ifdef MHD
                              const real Emag            = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg      );
                              const real Emag_IntT       = ( MagIntTime )
                                                         ? MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, MagSg_IntT )
                                                         : NULL_REAL;
#                             else
                              const real Emag            = NULL_REAL;
                              const real Emag_IntT       = NULL_REAL;
#                             endif
                              const real Eint = ( FluIntTime )
                                              ?   FluWeighting     *Hydro_Con2Eint( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    CheckMinEint_No, NULL_REAL, Emag )
                                                + FluWeighting_IntT*Hydro_Con2Eint( FluidPtr_IntT[DENS][k][j][i],
                                                                                    FluidPtr_IntT[MOMX][k][j][i],
                                                                                    FluidPtr_IntT[MOMY][k][j][i],
                                                                                    FluidPtr_IntT[MOMZ][k][j][i],
                                                                                    FluidPtr_IntT[ENGY][k][j][i],
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    CheckMinEint_No, NULL_REAL, Emag_IntT )
                                              :                     Hydro_Con2Eint( FluidPtr     [DENS][k][j][i],
                                                                                    FluidPtr     [MOMX][k][j][i],
                                                                                    FluidPtr     [MOMY][k][j][i],
                                                                                    FluidPtr     [MOMZ][k][j][i],
                                                                                    FluidPtr     [ENGY][k][j][i],
                                                                                    EoS_GuessHTilde_CPUPtr,
                                                                                    EoS_HTilde2Temp_CPUPtr,
                                                                                    CheckMinEint_No, NULL_REAL, Emag );
#                             endif // #ifdef DUAL_ENERGY ... else

                              OMP_Data  [p][TID][bin] += Eint*Weight;
                              OMP_Weight[p][TID][bin] += Weight;
                              OMP_NCell [p][TID][bin] ++;
                           }
                           break;
#                          endif // HYDRO

                           default:
                              Aux_Error( ERROR_INFO, "unsupported field (%ld) !!\n", TVarBitIdx[p] );
                              exit( 1 );
                        } // switch ( TVarBitIdx[p] )
                     } // if ( TFluIntIdx[p] != IdxUndef ) ... else ...
                  } // for (int p=0; p<NProf; p++)
               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=lv_min; lv<=lv_max; lv++)

#     if ( MODEL == HYDRO )
      delete [] Passive;         Passive      = NULL;
      delete [] Passive_IntT;    Passive_IntT = NULL;
#     endif

   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int p=0; p<NProf; p++)
   {
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         Prof[p]->Data  [b]  = OMP_Data  [p][0][b];
         Prof[p]->Weight[b]  = OMP_Weight[p][0][b];
         Prof[p]->NCell [b]  = OMP_NCell [p][0][b];
      }

      for (int t=1; t<NT; t++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         Prof[p]->Data  [b] += OMP_Data  [p][t][b];
         Prof[p]->Weight[b] += OMP_Weight[p][t][b];
         Prof[p]->NCell [b] += OMP_NCell [p][t][b];
      }
   }

// free per-thread arrays
   Aux_DeallocateArray3D( OMP_Data );
   Aux_DeallocateArray3D( OMP_Weight );
   Aux_DeallocateArray3D( OMP_NCell );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   for (int p=0; p<NProf; p++)
   {
      if ( MPI_Rank == 0 )
      {
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->Data,   Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->Weight, Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->NCell , Prof[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }

      else
      {
         MPI_Reduce( Prof[p]->Data,   NULL,            Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Prof[p]->Weight, NULL,            Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Prof[p]->NCell,  NULL,            Prof[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }
   }
#  endif


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   {
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
//       skip empty bins since both their data and weight are zero
         if ( Prof[p]->NCell[b] > 0L )    Prof[p]->Data[b] /= Prof[p]->Weight[b];
      }
   }


// broadcast data to all ranks
   for (int p=0; p<NProf; p++)
   {
      MPI_Bcast( Prof[p]->Data,   Prof[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Prof[p]->Weight, Prof[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Prof[p]->NCell,  Prof[p]->NBin, MPI_LONG,   0, MPI_COMM_WORLD );
   }


// remove the empty bins
// --> all ranks do the same work so that no data broadcast is required
   if ( RemoveEmpty )
   {
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         if ( Prof[0]->NCell[b] != 0L )   continue;

//       remove consecutive empty bins at the same time for better performance
         int b_up;
         for (b_up=b+1; b_up<Prof[0]->NBin; b_up++)
            if ( Prof[0]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (b_up=b+stride; b_up<Prof[0]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            for (int p=0; p<NProf; p++)
            {
               Prof[p]->Radius[b_up_ms] = Prof[p]->Radius[b_up];
               Prof[p]->Data  [b_up_ms] = Prof[p]->Data  [b_up];
               Prof[p]->Weight[b_up_ms] = Prof[p]->Weight[b_up];
               Prof[p]->NCell [b_up_ms] = Prof[p]->NCell [b_up];
            }
         }

//       reset the total number of bins
         for (int p=0; p<NProf; p++)
            Prof[p]->NBin -= stride;
      } // for (int b=0; b<Prof->NBin; b++)

//    update the maximum radius since the last bin may have been removed
      for (int p=0; p<NProf; p++)
      {
         const int LastBin = Prof[p]->NBin-1;

         Prof[p]->MaxRadius = ( LogBin ) ? Prof[p]->Radius[LastBin]*sqrt( LogBinRatio )
                                         : Prof[p]->Radius[LastBin] + 0.5*dr_min;
      }
   } // if ( RemoveEmpty )

} // FUNCTION : Aux_ComputeProfile
