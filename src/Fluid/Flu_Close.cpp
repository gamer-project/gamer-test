#include "GAMER.h"
#include "CUFLU.h"


// status of the fluid solver used by AUTO_REDUCE_DT (declared in Flu_AdvanceDt.cpp)
extern int FluStatus_ThisRank;

// whether or not to continue applying AUTO_REDUCE_DT (decalred in Flu_AdvanceDt.cpp)
extern bool AutoReduceDt_Continue;


static void StoreFlux( const int lv, const real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                       const int NPG, const int *PID0_List, const real dt );
static void CorrectFlux( const int SonLv, const real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                         const int NPG, const int *PID0_List, const real dt );
#if ( MODEL == HYDRO )
static bool Unphysical( const real Fluid[], const int CheckMode, const real Emag );
static void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                               const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                               real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                               char h_DE_Array_F_Out[][ CUBE(PS2) ],
                               real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                               const real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                               const real h_Mag_Array_F_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                               const real dt );
#ifdef MHD
void StoreElectric( const int lv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                    const int NPG, const int *PID0_List, const real dt );
void CorrectElectric( const int SonLv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                      const int NPG, const int *PID0_List, const real dt );
#endif
#endif // #if ( MODEL == HYDRO )
extern void Hydro_RiemannSolver_Roe ( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                      const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                                      const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
extern void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                      const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                                      const EoS_H2TEM_t EoS_HTilde2Temp, const EoS_TEM2C_t EoS_Temper2CSqr,
                                      const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                                      const real* const EoS_Table[EOS_NTABLE_MAX] );
extern void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                      const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                                      const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Close
// Description :  1. Save the fluxes across the coarse-fine boundaries at level "lv"
//                2. Correct the fluxes across the coarse-fine boundaries at level "lv-1"
//                3. Copy the data from the "h_Flu_Array_F_Out" and "h_DE_Array_F_Out" arrays to the "amr->patch" pointers
//                4. Get the minimum time-step information of the fluid solver
//
// Parameter   :  lv                : Target refinement level
//                SaveSg_Flu        : Sandglass to store the updated fluid data
//                SaveSg_Mag        : Sandglass to store the updated B field (for MHD only)
//                h_Flux_Array      : Host array storing the updated flux data
//                h_Ele_Array       : Host array storing the updated E field (for MHD only)
//                h_Flu_Array_F_Out : Host array storing the updated fluid data
//                h_Mag_Array_F_Out : Host array storing the updated B field (for MHD only)
//                h_DE_Array_F_Out  : Host array storing the dual-energy status (for DUAL_ENERGY only)
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indices with LocalID==0 to be udpated
//                h_Flu_Array_F_In  : Host array storing the input fluid variables
//                h_Mag_Array_F_In  : Host array storing the input B field (for MHD only)
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void Flu_Close( const int lv, const int SaveSg_Flu, const int SaveSg_Mag,
                real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                real h_Mag_Array_F_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                char h_DE_Array_F_Out[][ CUBE(PS2) ],
                const int NPG, const int *PID0_List,
                const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                const real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                const double dt )
{

// try to correct the unphysical results in h_Flu_Array_F_Out (e.g., negative density)
// --> must be done BEFORE invoking both StoreFlux() and CorrectFlux() since CorrectUnphysical() might modify the flux array
#  if ( MODEL == HYDRO && !defined SRHD )
   CorrectUnphysical( lv, NPG, PID0_List, h_Flu_Array_F_In, h_Flu_Array_F_Out, h_DE_Array_F_Out, h_Flux_Array,
                      h_Mag_Array_F_In, h_Mag_Array_F_Out, dt );
#  endif


// operations related to flux fix-up
   if ( OPT__FIXUP_FLUX )
   {
//    save the fluxes on the coarse-fine boundaries at level "lv"
      StoreFlux( lv, h_Flux_Array, NPG, PID0_List, dt );

//    correct the fluxes on the coarse-fine boundaries at level "lv-1"
      CorrectFlux( lv, h_Flux_Array, NPG, PID0_List, dt );
   }


// operations related to electric field fix-up
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )
   {
//    save the E field on the coarse-fine boundaries at level "lv"
      if ( lv != TOP_LEVEL )  StoreElectric( lv, h_Ele_Array, NPG, PID0_List, dt );

//    correct the E field on the coarse-fine boundaries at level "lv-1"
      CorrectElectric( lv, h_Ele_Array, NPG, PID0_List, dt );
   }
#  endif


// copy the updated data from output arrays to the corresponding patch pointers
#  if ( FLU_NOUT != NCOMP_TOTAL )
#     error : ERROR : FLU_NOUT != NCOMP_TOTAL (one must specify how to copy data from h_Flu_Array_F_Out to fluid) !!
#  endif

#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID     = PID0 + LocalID;
         const int Table_x = TABLE_02( LocalID, 'x', 0, PATCH_SIZE );
         const int Table_y = TABLE_02( LocalID, 'y', 0, PATCH_SIZE );
         const int Table_z = TABLE_02( LocalID, 'z', 0, PATCH_SIZE );

         int I, J, K, KJI;

//       fluid variables
         for (int v=0; v<FLU_NOUT; v++)      {
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = IDX321( I, J, K, PS2, PS2 );

            amr->patch[SaveSg_Flu][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_F_Out[TID][v][KJI];

         }}}}

//       dual-energy status
#        ifdef DUAL_ENERGY
         for (int k=0; k<PATCH_SIZE; k++)    {  K = Table_z + k;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = Table_y + j;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = Table_x + i;

            KJI = IDX321( I, J, K, PS2, PS2 );

//          de_status is always stored in Sg=0
            amr->patch[0][lv][PID]->de_status[k][j][i] = h_DE_Array_F_Out[TID][KJI];

         }}}
#        endif

//       magnetic field
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG; v++)
         {
            int ijk_end[3], idx=0;

            for (int d=0; d<3; d++)    ijk_end[d] = ( d == v ) ? PS1+1 : PS1;

            for (int k=0; k<ijk_end[2]; k++)    {  K = Table_z + k;
            for (int j=0; j<ijk_end[1]; j++)    {  J = Table_y + j;
            for (int i=0; i<ijk_end[0]; i++)    {  I = Table_x + i;

               KJI = IDX321( I, J, K, ijk_end[0]+PS1, ijk_end[1]+PS1 );

               amr->patch[SaveSg_Mag][lv][PID]->magnetic[v][ idx ++ ] = h_Mag_Array_F_Out[TID][v][KJI];

            }}}
         } // for (int v=0; v<NCOMP_MAG; v++)
#        endif // #ifdef MHD

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Flu_Close



//-------------------------------------------------------------------------------------------------------
// Function    :  StoreFlux
// Description :  Save the coarse-grid fluxes across the coarse-fine boundaries for patches at level "lv"
//                --> to be fixed later by CorrectFlux()
//
// Parameter   :  lv           : Target refinement level
//                h_Flux_Array : Host array storing the updated flux data
//                NPG          : Number of patch groups to be evaluated
//                PID0_List    : List recording the patch indices with LocalID==0 to be udpated
//                dt           : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void StoreFlux( const int lv, const real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )   Aux_Error( ERROR_INFO, "amr->WithFlux is off !!\n" );


// nothing to do on the highest level
   if ( lv == TOP_LEVEL )  return;


#  pragma omp parallel for schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;

         for (int s=0; s<6; s++)
         {
//          for bitwise reproducibility, store the fluxes to be corrected in flux_bitrep[]
#           ifdef BIT_REP_FLUX
            real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux_bitrep[s];
#           else
            real (*FluxPtr)[PS1][PS1] = amr->patch[0][lv][PID]->flux[s];
#           endif

            if ( FluxPtr != NULL )
            {
               int face_idx, disp_m, disp_n;

               switch ( s )
               {
                  case 0:  case 1:
                     face_idx = TABLE_02( LocalID, 'x', 0, 1 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'y', 0, PS1 );
                     break;

                  case 2:  case 3:
                     face_idx = TABLE_02( LocalID, 'y', 1, 2 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  case 4:  case 5:
                     face_idx = TABLE_02( LocalID, 'z', 2, 3 ) + s;
                     disp_m   = TABLE_02( LocalID, 'y', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
               }


               for (int v=0; v<NFLUX_TOTAL; v++)   {
               for (int m=0; m<PS1; m++)           {  const int mm = m + disp_m;
               for (int n=0; n<PS1; n++)           {  const int nn = n + disp_n;

                  FluxPtr[v][m][n] = dt*h_Flux_Array[TID][face_idx][v][ mm*PS2 + nn ];

               }}}

            } // if ( FluxPtr != NULL )
         } // for (int s=0; s<6; s++)
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : StoreFlux



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectFlux
// Description :  Use the fluxes across the coarse-fine boundaries at level SonLv to correct the fluxes
//                at level FaLv
//
// Parameter   :  SonLv        : Target refinement level
//                h_Flux_Array : Array storing the updated flux data
//                NPG          : Number of patch groups to be evaluated
//                PID0_List    : List recording the patch indices with LocalID==0 to be udpated
//                dt           : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectFlux( const int SonLv, const real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                  const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithFlux )   Aux_Error( ERROR_INFO, "amr->WithFlux is off !!\n" );


// nothing to do on the root level
   if ( SonLv == 0 )    return;


   const int  FaLv         = SonLv - 1;
   const int  Mapping  [6] = { 0, 2, 3, 5, 6, 8 };
   const int  MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
   const real dt_4         = (real)0.25*dt;

#  pragma omp parallel
   {
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         const int PID0  = PID0_List[TID];
         const int FaPID = amr->patch[0][SonLv][PID0]->father;

#        ifdef GAMER_DEBUG
         if ( FaPID < 0 )
            Aux_Error( ERROR_INFO, "FaPID = %d < 0 (SonLv %d, PID0 %d) !!\n", FaPID, SonLv, PID0 );
#        endif

         for (int s=0; s<6; s++)
         {
            const int FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#           ifdef GAMER_DEBUG
            if ( FaSibPID == -1 )
               Aux_Error( ERROR_INFO, "FaSibPID == -1 (FaLv %d, FaPID %d, s %d, PID0 %d) !!\n", FaLv, FaPID, s, PID0 );
#           endif

//          skip patches adjacent to non-periodic boundaries
            if ( FaSibPID < -1 )    continue;

//          for AUTO_REDUCE_DT, store the updated fluxes in the temporary array flux_tmp[] since
//          we may need to abandon them if the fluid solver fails
            real (*FluxPtr)[PS1][PS1] = ( AUTO_REDUCE_DT ) ? amr->patch[0][FaLv][FaSibPID]->flux_tmp[ MirrorSib[s] ] :
                                                             amr->patch[0][FaLv][FaSibPID]->flux    [ MirrorSib[s] ];

//          skip patches not adjacent to coarse-fine boundaries
            if ( FluxPtr == NULL )  continue;

//          store fluxes
            for (int v=0; v<NFLUX_TOTAL; v++)
            for (int m=0; m<PS2; m++)  { const int mm = m/2;
            for (int n=0; n<PS2; n++)  { const int nn = n/2;

               FluxPtr[v][mm][nn] -= dt_4*h_Flux_Array[TID][ Mapping[s] ][v][ m*PS2 + n ];

            }}
         } // for (int s=0; s<6; s++)
      } // for (int TID=0; TID<NPG; TID++)

   } // OpenMP parallel region

} // FUNCTION : CorrectFlux



#if ( MODEL == HYDRO &&  !defined SRHD )
//-------------------------------------------------------------------------------------------------------
// Function    :  Unphysical
// Description :  Check whether the input variables are unphysical
//
// Note        :  1. One can put arbitrary criteria here. Cells violating the conditions will be recalculated
//                   in CorrectUnphysical()
//                2. Currently only used for HYDRO to check whether the input mass density, internal energy/pressure,
//                   or total energy is smaller than the given thresholds
//                   --> It also checks if any variable is -inf, +inf, or nan
//                   --> It does NOT check if passive scalars are negative
//                       --> We already apply a floor value in Hydro_Shared_FullStepUpdate()
//                3. When enabling the dual-energy formalism (with DE_ENPY), CheckMode=1 checks pressure
//                   instead of internal energy
//
// Parameter   :  Fluid     : Input fluid variable array with size FLU_NOUT
//                CheckMode : (0/1) --> check (total energy/internal energy)
//                Emag      : Magnetic energy (for MHD only)
//
// Return      :  true/false <--> input Fluid[] is unphysical/physical
//-------------------------------------------------------------------------------------------------------
bool Unphysical( const real Fluid[], const int CheckMode, const real Emag )
{

   const int  CheckMinEtot = 0;
   const int  CheckMinEint = 1;
   const bool NoFloor      = false;


// if any checks below fail, return true
// =================================================
// note that since MIN_DENS and MIN_PRES are declared as double, they must be converted to **real** before the comparison
// --> otherwise LHS in the comparison will be converted from real to double, which is inconsistent with the assignment
//     (e.g., "Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS" )
   if (  !Aux_IsFinite(Fluid[DENS])  ||  !Aux_IsFinite(Fluid[MOMX])  ||  !Aux_IsFinite(Fluid[MOMY])  ||
         !Aux_IsFinite(Fluid[MOMZ])  ||  !Aux_IsFinite(Fluid[ENGY])  ||
         Fluid[DENS] < (real)MIN_DENS  )
      return true;

#  ifndef BAROTROPIC_EOS
   if ( CheckMode == CheckMinEtot  &&  Fluid[ENGY] < (real)MIN_EINT )
      return true;

   if ( CheckMode == CheckMinEint  &&
         (
//          when adopting the dual-energy formalism, do NOT calculate pressure from "Etot-Ekin" since it would suffer
//          from large round-off errors
//          --> currently we use TINY_NUMBER as entropy floor and hence here we use 2.0*TINY_NUMBER to
//              validate entropy
//          --> in general, MIN_PRES > 0.0 should be sufficient for detecting unphysical entropy
//          --> however, the additional check "Fluid[ENPY] < (real)2.0*TINY_NUMBER" is necessary when MIN_PRES == 0.0
#           if   ( DUAL_ENERGY == DE_ENPY )
            Hydro_DensEntropy2Pres( Fluid[DENS], Fluid[ENPY], EoS_AuxArray_Flt[1], NoFloor, NULL_REAL ) < (real)MIN_PRES  ||
            Fluid[ENPY] < (real)2.0*TINY_NUMBER

#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!

#           else // without DUAL_ENERGY
            Hydro_Con2Eint( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                            NULL, NULL, NoFloor, NULL_REAL, Emag ) < (real)MIN_EINT

#           endif // DUAL_ENERGY
         )
      )
      return true;
#  endif // #ifndef BAROTROPIC_EOS


// if all checks above pass, return false
// =================================================
   return false;

} // FUNCTION : Unphysical



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectUnphysical
// Description :  Check if any cell in the output array of Fluid solver contains unphysical results
//                --> For example, negative density
//
// Note        :  1. Define unphysical values in Unphysical()
//                2. Currently only used for HYDRO to check whether the input mass density, internal energy/pressure,
//                   or total energy is smaller than the minimum allowed values (i.e., MIN_DENS/MIN_PRES/MIN_EINT)
//                   --> It also checks if any variable is -inf, +inf, and nan
//                   --> But one can define arbitrary criteria in Unphysical() to trigger the correction
//                3. Procedure:
//                   if ( found_unphysical )
//                   {
//                      if ( OPT__1ST_FLUX_CORR )
//                         Try to correct unphysical variables using 1st-order fluxes
//                      else
//                         Do nothing
//
//                      if ( ! AUTO_REDUCE_DT )
//                         Apply floors
//
//                      if ( still_found_unphysical )
//                         if ( AUTO_REDUCE_DT )
//                            Invoke the fluid solver again on the same level but with a smaller dt
//                         else
//                            Print debug messages and abort
//                      else
//                         Store the corrected results in the output array
//                   }
//
// Parameter   :  lv                : Target refinement level
//                NPG               : Number of patch groups to be evaluated
//                PID0_List         : List recording the patch indices with LocalID==0 to be udpated (for debug only)
//                h_Flu_Array_F_In  : Input fluid array
//                h_Flu_Array_F_Out : Output fluid array
//                h_DE_Array_F_Out  : Output dual-energy status array
//                h_Flux_Array      : Output array storing the updated flux data
//                h_Mag_Array_F_In  : Input B field array
//                h_Mag_Array_F_Out : Output B field array
//                dt                : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectUnphysical( const int lv, const int NPG, const int *PID0_List,
                        const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                        real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                        char h_DE_Array_F_Out[][ CUBE(PS2) ],
                        real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                        const real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                        const real h_Mag_Array_F_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                        const real dt )
{

   const real dh               = (real)amr->dh[lv];
   const real dt_dh            = dt/dh;
   const int  didx[3]          = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   const int  CheckMinEtot     = 0;
   const int  CheckMinEint     = 1;
   const int  LocalID[2][2][2] = { 0, 1, 2, 4, 3, 6, 5, 7 };

// variables shared by all OpenMP threads
   long NCorrThisTime = 0;
   bool CorrectUnphy  = GAMER_SUCCESS;


// OpenMP parallel region
#  pragma omp parallel
   {

// variables private to each OpenMP thread
   real VarL[3][NCOMP_TOTAL], VarC[NCOMP_TOTAL], VarR[3][NCOMP_TOTAL], FluxL[3][NCOMP_TOTAL], FluxR[3][NCOMP_TOTAL];
   real dF[3][NCOMP_TOTAL], Out[NCOMP_TOTAL], Update[NCOMP_TOTAL];
   real FluxL_1D[NCOMP_TOTAL], FluxR_1D[NCOMP_TOTAL];
   int  ijk_out[3];

// variables for OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D
   const int  Corr1D_NBuf          = 1;
   const int  Corr1D_NCell         = 2*Corr1D_NBuf + 1;
// [3][3] = [dimension][i/j/k indices]
// we apply the directionally splitting correction in the order x->y->z
   const int  Corr1D_idx_min[3][3] = { { Corr1D_NBuf,              0,              0 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,              0 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,    Corr1D_NBuf } };
   const int  Corr1D_idx_max[3][3] = { { Corr1D_NBuf, Corr1D_NCell-1, Corr1D_NCell-1 },
                                       { Corr1D_NBuf,    Corr1D_NBuf, Corr1D_NCell-1 },
                                       { Corr1D_NBuf,    Corr1D_NBuf,    Corr1D_NBuf } };
   const int  Corr1D_didx1[3]      = { NCOMP_TOTAL, Corr1D_NCell*NCOMP_TOTAL, SQR(Corr1D_NCell)*NCOMP_TOTAL };
#  if ( DUAL_ENERGY == DE_ENPY )
   const bool CorrPres_Yes         = true;
   const bool CorrPres_No          = false;
#  endif

   int   Corr1D_didx2[3];
   real *Corr1D_InOut_PtrL=NULL, *Corr1D_InOut_PtrC=NULL, *Corr1D_InOut_PtrR=NULL;
   real (*Corr1D_InOut)[Corr1D_NCell][Corr1D_NCell][NCOMP_TOTAL]
         = ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D ) ? new real [Corr1D_NCell][Corr1D_NCell][Corr1D_NCell][NCOMP_TOTAL]
                                                          : NULL;

#  pragma omp for reduction( +:NCorrThisTime ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      for (ijk_out[2]=0; ijk_out[2]<PS2; ijk_out[2]++)
      for (ijk_out[1]=0; ijk_out[1]<PS2; ijk_out[1]++)
      for (ijk_out[0]=0; ijk_out[0]<PS2; ijk_out[0]++)
      {
         const int idx_out = (ijk_out[2]*PS2 + ijk_out[1])*PS2 + ijk_out[0];

//       check if the updated values are unphysical
         for (int v=0; v<NCOMP_TOTAL; v++)   Out[v] = h_Flu_Array_F_Out[TID][v][idx_out];

//       compute the magnetic energy
#        ifdef MHD
         const real Emag_Out = MHD_GetCellCenteredBEnergy( h_Mag_Array_F_Out[TID][MAGX],
                                                           h_Mag_Array_F_Out[TID][MAGY],
                                                           h_Mag_Array_F_Out[TID][MAGZ],
                                                           PS2, PS2, PS2, ijk_out[0], ijk_out[1], ijk_out[2] );
#        else
         const real Emag_Out = NULL_REAL;
#        endif

         if ( Unphysical(Out, CheckMinEint, Emag_Out) )
         {
            const int idx_in = ( (ijk_out[2]+FLU_GHOST_SIZE)*FLU_NXT + (ijk_out[1]+FLU_GHOST_SIZE) )*FLU_NXT + (ijk_out[0]+FLU_GHOST_SIZE);

//          try to correct unphysical results using 1st-order fluxes (which should be more diffusive)
//          ========================================================================================
//          (A) first apply the directionally **unsplitting** correction
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
            {
//             collect nearby input coserved variables
               for (int v=0; v<NCOMP_TOTAL; v++)
               {
//                here we have assumed that the fluid solver does NOT modify the input fluid array
//                --> not applicable to RTVD
                  VarC[v] = h_Flu_Array_F_In[TID][v][idx_in];

                  for (int d=0; d<3; d++)
                  {
                     VarL[d][v] = h_Flu_Array_F_In[TID][v][ idx_in - didx[d] ];
                     VarR[d][v] = h_Flu_Array_F_In[TID][v][ idx_in + didx[d] ];
                  }
               }

//             invoke Riemann solver to calculate the fluxes
//             (note that the recalculated flux does NOT include gravity even for UNSPLIT_GRAVITY --> reduce to 1st-order accuracy)
               switch ( OPT__1ST_FLUX_CORR_SCHEME )
               {
                  case RSOLVER_1ST_ROE:
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_Roe ( d, FluxL[d], VarL[d], VarC,    MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        Hydro_RiemannSolver_Roe ( d, FluxR[d], VarC,    VarR[d], MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     }
                     break;

#                 ifndef MHD
                  case RSOLVER_1ST_HLLC:
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_HLLC( d, FluxL[d], VarL[d], VarC,    MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temper2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        Hydro_RiemannSolver_HLLC( d, FluxR[d], VarC,    VarR[d], MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temper2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     }
                     break;
#                 endif

                  case RSOLVER_1ST_HLLE:
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_HLLE( d, FluxL[d], VarL[d], VarC,    MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        Hydro_RiemannSolver_HLLE( d, FluxR[d], VarC,    VarR[d], MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     }
                     break;

#                 ifdef MHD
                  case RSOLVER_1ST_HLLD:
                     Aux_Error( ERROR_INFO, "RSOLVER_1ST_HLLD in MHD is NOT supported yet !!\n" );
                     /*
                     for (int d=0; d<3; d++)
                     {
                        Hydro_RiemannSolver_HLLD( d, FluxL[d], VarL[d], VarC,    MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                        Hydro_RiemannSolver_HLLD( d, FluxR[d], VarC,    VarR[d], MIN_DENS, MIN_PRES,
                                                  EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                     }
                     */
                     break;
#                 endif

                  default:
                     Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__1ST_FLUX_CORR_SCHEME );
               } // switch ( OPT__1ST_FLUX_CORR_SCHEME )

//             recalculate the first-order solution for a full time-step
               for (int d=0; d<3; d++)
               for (int v=0; v<NCOMP_TOTAL; v++)   dF[d][v] = FluxR[d][v] - FluxL[d][v];

               for (int v=0; v<NCOMP_TOTAL; v++)
                  Update[v] = h_Flu_Array_F_In[TID][v][idx_in] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );

            } // if ( OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )


//          ========================================================================================
//          (B) apply the directionally **splitting** correction if the unsplitting correction failed
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D )
            {
//             apply the dual-energy formalism to correct the internal energy
//             --> gravity solver may update the internal energy and dual-energy variable again when
//                 UNSPLIT_GRAVITY is adopted
//             --> if the corrected internal energy (pressure) passes the test, we don't need to apply the
//                 directionally **splitting** correction
//             --> we do NOT apply the minimum pressure check in Hydro_DualEnergyFix() here
//                 --> otherwise the pressure floor might disable the 1st-order-flux correction
#              ifdef DUAL_ENERGY
               Hydro_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                                    h_DE_Array_F_Out[TID][idx_out], EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2],
                                    CorrPres_No, NULL_REAL, DUAL_ENERGY_SWITCH, Emag_Out );
#              endif

               if ( Unphysical(Update, CheckMinEint, Emag_Out) )
               {
//                collect nearby input coserved variables
                  for (int k=0; k<Corr1D_NCell; k++)  { Corr1D_didx2[2] = (k-Corr1D_NBuf)*didx[2];
                  for (int j=0; j<Corr1D_NCell; j++)  { Corr1D_didx2[1] = (j-Corr1D_NBuf)*didx[1];
                  for (int i=0; i<Corr1D_NCell; i++)  { Corr1D_didx2[0] = (i-Corr1D_NBuf)*didx[0];

//                   here we have assumed that the fluid solver does NOT modify the input fluid array
//                   --> not applicable to RTVD
                     for (int v=0; v<NCOMP_TOTAL; v++)
                        Corr1D_InOut[k][j][i][v] = h_Flu_Array_F_In[TID][v][ idx_in + Corr1D_didx2[0] + Corr1D_didx2[1] + Corr1D_didx2[2] ];
                  }}}

//                apply the 1st-order-flux correction along different spatial directions **successively**
                  for (int d=0; d<3; d++)
                  {
                     for (int k=Corr1D_idx_min[d][2]; k<=Corr1D_idx_max[d][2]; k++)
                     for (int j=Corr1D_idx_min[d][1]; j<=Corr1D_idx_max[d][1]; j++)
                     for (int i=Corr1D_idx_min[d][0]; i<=Corr1D_idx_max[d][0]; i++)
                     {
//                      get the pointers to the left and right states for the Riemann solvers
                        Corr1D_InOut_PtrC = Corr1D_InOut[k][j][i];
                        Corr1D_InOut_PtrL = Corr1D_InOut_PtrC - Corr1D_didx1[d];
                        Corr1D_InOut_PtrR = Corr1D_InOut_PtrC + Corr1D_didx1[d];

//                      invoke Riemann solver to calculate the fluxes
//                      (note that the recalculated flux does NOT include gravity even for UNSPLIT_GRAVITY --> reduce to 1st-order accuracy)
                        switch ( OPT__1ST_FLUX_CORR_SCHEME )
                        {
                           case RSOLVER_1ST_ROE:
                              Hydro_RiemannSolver_Roe ( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                              Hydro_RiemannSolver_Roe ( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                           break;

#                          ifndef MHD
                           case RSOLVER_1ST_HLLC:
                              Hydro_RiemannSolver_HLLC( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temper2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                              Hydro_RiemannSolver_HLLC( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temper2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                           break;
#                          endif

                           case RSOLVER_1ST_HLLE:
                              Hydro_RiemannSolver_HLLE( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                              Hydro_RiemannSolver_HLLE( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                           break;

#                          ifdef MHD
                           case RSOLVER_1ST_HLLD:
                              Aux_Error( ERROR_INFO, "RSOLVER_1ST_HLLD in MHD is NOT supported yet !!\n" );
                              /*
                              Hydro_RiemannSolver_HLLD( d, FluxL_1D, Corr1D_InOut_PtrL, Corr1D_InOut_PtrC, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                              Hydro_RiemannSolver_HLLD( d, FluxR_1D, Corr1D_InOut_PtrC, Corr1D_InOut_PtrR, MIN_DENS, MIN_PRES,
                                                        EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                              */
                           break;
#                          endif

                           default:
                              Aux_Error( ERROR_INFO, "unnsupported Riemann solver (%d) !!\n", OPT__1ST_FLUX_CORR_SCHEME );
                        }

//                      recalculate the first-order solution for a full time-step
                        for (int v=0; v<NCOMP_TOTAL; v++)
                           Corr1D_InOut[k][j][i][v] -= dt_dh*( FluxR_1D[v] - FluxL_1D[v] );

//                      store the 1st-order fluxes used for updating the central cell
                        if ( i == Corr1D_NBuf  &&  j == Corr1D_NBuf  &&  k == Corr1D_NBuf )
                        {
                           for (int v=0; v<NCOMP_TOTAL; v++)
                           {
                              FluxL[d][v] = FluxL_1D[v];
                              FluxR[d][v] = FluxR_1D[v];
                           }
                        }
                     } // i,j,k
                  } // for (int d=0; d<3; d++)

//                store the corrected results in Update[]
                  for (int v=0; v<NCOMP_TOTAL; v++)
                     Update[v] = Corr1D_InOut[Corr1D_NBuf][Corr1D_NBuf][Corr1D_NBuf][v];

               } // if ( Unphysical(Update, CheckMinEint, Emag_Out) )
            } // if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D )


//          ========================================================================================
//          (C) if OPT__1ST_FLUX_CORR is off, just copy data to Update[] for checking unphysical results
//              in the next step
//          ========================================================================================
            if ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_NONE )
            {
               for (int v=0; v<NCOMP_TOTAL; v++)   Update[v] = Out[v];
            }


//          ensure positive density
//          --> apply it only when AutoReduceDt_Continue is false
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this density floor
//          --> note that MIN_DENS is declared as double and must be converted to **real** before the comparison
//              --> to be consistent with the check in Unphysical()
//          --> do NOT check the minimum internal energy here since we want to apply the dual-energy correction first
            if ( ! AutoReduceDt_Continue  &&  OPT__LAST_RESORT_FLOOR )
               Update[DENS] = FMAX( Update[DENS], (real)MIN_DENS );


//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Update[v] = FMAX( Update[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               Hydro_NormalizePassive( Update[DENS], Update+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif


//          apply the dual-energy formalism to correct the internal energy again
//          --> gravity solver may update the internal energy and dual-energy variable again when
//              UNSPLIT_GRAVITY is adopted
//          --> this might be redundant when OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D
//              --> but it ensures the consistency between all fluid variables since we apply density floor AFTER
//                  the 1st-order-flux correction
//          --> we apply the minimum pressure check in Hydro_DualEnergyFix() here only when AutoReduceDt_Continue is false
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this pressure floor
#           ifdef DUAL_ENERGY
            Hydro_DualEnergyFix( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY], Update[ENPY],
                                 h_DE_Array_F_Out[TID][idx_out], EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2],
                                 (!AutoReduceDt_Continue && OPT__LAST_RESORT_FLOOR) ? CorrPres_Yes : CorrPres_No,
                                 MIN_PRES, DUAL_ENERGY_SWITCH, Emag_Out );

//          apply internal energy floor if dual-energy formalism is not adopted
//          --> apply it only when AutoReduceDt_Continue is false
//              --> otherwise AUTO_REDUCE_DT may not be triggered due to this internal energy floor
#           else
            if ( ! AutoReduceDt_Continue  &&  OPT__LAST_RESORT_FLOOR )
               Update[ENGY] = Hydro_CheckMinEintInEngy( Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                                        MIN_EINT, Emag_Out );
#           endif


//          check if the newly updated values are still unphysical
//          --> note that, when AutoReduceDt_Continue is false, we check Etot instead of Eint since even after calling
//              Hydro_CheckMinEintInEngy() we may still have Eint < MIN_EINT due to round-off errors (especially when Eint << Ekin)
//              --> it will not crash the code since we always apply MIN_EINT/MIN_PRES when calculating Eint/pressure
//          --> when AutoReduceDt_Continue is true, we still check Eint instead of Etot
            if ( Unphysical(Update, (AutoReduceDt_Continue)?CheckMinEint:CheckMinEtot, Emag_Out) )
            {
//             set CorrectUnphy = GAMER_FAILED if any cells fail
//             --> use critical directive to avoid thread racing (may not be necessary here?)
#              pragma omp critical
               CorrectUnphy = GAMER_FAILED;


//             output the debug information (only if AutoReduceDt_Continue is false)
               if ( ! AutoReduceDt_Continue )
               {
                  const int  PID_Failed      = PID0_List[TID] + LocalID[ijk_out[2]/PS1][ijk_out[1]/PS1][ijk_out[0]/PS1];
                  const bool CheckMinEint_No = false;
                  real In[NCOMP_TOTAL], tmp[NCOMP_TOTAL];

                  char FileName[100];
                  sprintf( FileName, "FailedPatchGroup_r%03d_lv%02d_PID0-%05d", MPI_Rank, lv, PID0_List[TID] );

//                use "a" instead of "w" since there may be more than one failed cell in a given patch group
                  FILE *File = fopen( FileName, "a" );

                  for (int v=0; v<NCOMP_TOTAL; v++)   In[v] = h_Flu_Array_F_In[TID][v][idx_in];

#                 ifdef MHD
                  const real Emag_In     = MHD_GetCellCenteredBEnergy( h_Mag_Array_F_In[TID][MAGX],
                                                                       h_Mag_Array_F_In[TID][MAGY],
                                                                       h_Mag_Array_F_In[TID][MAGZ],
                                                                       FLU_NXT, FLU_NXT, FLU_NXT,
                                                                       ijk_out[0]+FLU_GHOST_SIZE,
                                                                       ijk_out[1]+FLU_GHOST_SIZE,
                                                                       ijk_out[2]+FLU_GHOST_SIZE );
                  const real Emag_Update = Emag_Out;
#                 else
                  const real Emag_In     = NULL_REAL;
                  const real Emag_Update = NULL_REAL;
#                 endif

//                output information about the failed cell
                  fprintf( File, "PID                              = %5d\n", PID_Failed );
                  fprintf( File, "(i,j,k) in the patch             = (%2d,%2d,%2d)\n", ijk_out[0]%PS1, ijk_out[1]%PS1, ijk_out[2]%PS1 );
                  fprintf( File, "(i,j,k) in the input fluid array = (%2d,%2d,%2d)\n", ijk_out[0], ijk_out[1], ijk_out[2] );
#                 ifdef DUAL_ENERGY
                  fprintf( File, "total energy density update      = %s\n",
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_ETOT     ) ? "Etot" :
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_DUAL     ) ? "Dual" :
                           ( h_DE_Array_F_Out[TID][idx_out] == DE_UPDATED_BY_MIN_PRES ) ? "MinPres" :
                                                                                          "Unknown" );
#                 endif
                  fprintf( File, "\n" );

                  fprintf( File, "               (%14s, %14s, %14s, %14s, %14s, %14s",
                           FieldLabel[DENS], FieldLabel[MOMX], FieldLabel[MOMY], FieldLabel[MOMZ], FieldLabel[ENGY], "Eint" );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14s", FieldLabel[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "input        = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                           Hydro_Con2Eint(In[DENS], In[MOMX], In[MOMY], In[MOMZ], In[ENGY],
                                          NULL, NULL, CheckMinEint_No, NULL_REAL, Emag_In) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", In[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "ouptut (old) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                           Hydro_Con2Eint(Out[DENS], Out[MOMX], Out[MOMY], Out[MOMZ], Out[ENGY],
                                          NULL, NULL, CheckMinEint_No, NULL_REAL, Emag_Out) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Out[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

                  fprintf( File, "output (new) = (%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e",
                           Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                           Hydro_Con2Eint(Update[DENS], Update[MOMX], Update[MOMY], Update[MOMZ], Update[ENGY],
                                          NULL, NULL, CheckMinEint_No, NULL_REAL, Emag_Update) );
#                 if ( DUAL_ENERGY == DE_ENPY )
                  fprintf( File, ", %14.7e", Update[ENPY] );
#                 endif
                  fprintf( File, ")\n" );

//                output all data in the input fluid array (including ghost zones)
                  fprintf( File, "\nFull input array including ghost zones\n" );
                  fprintf( File, "===============================================================================================\n" );
                  fprintf( File, "(%2s,%2s,%2s)", "i", "j", "k" );
                  for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %14s", FieldLabel[v] );
                  fprintf( File, " %14s", "Eint" );
                  fprintf( File, "\n" );

                  for (int k=0; k<FLU_NXT; k++)
                  for (int j=0; j<FLU_NXT; j++)
                  for (int i=0; i<FLU_NXT; i++)
                  {
                     fprintf( File, "(%2d,%2d,%2d)", i-FLU_GHOST_SIZE, j-FLU_GHOST_SIZE, k-FLU_GHOST_SIZE );

                     for (int v=0; v<NCOMP_TOTAL; v++)
                     {
                        tmp[v] = h_Flu_Array_F_In[TID][v][ ((k*FLU_NXT)+j)*FLU_NXT+i ];
                        fprintf( File, " %14.7e", tmp[v] );
                     }

#                    ifdef MHD
                     const real Emag_tmp = MHD_GetCellCenteredBEnergy( h_Mag_Array_F_In[TID][MAGX],
                                                                       h_Mag_Array_F_In[TID][MAGY],
                                                                       h_Mag_Array_F_In[TID][MAGZ],
                                                                       FLU_NXT, FLU_NXT, FLU_NXT, i, j, k );
#                    else
                     const real Emag_tmp = NULL_REAL;
#                    endif

                     fprintf( File, " %14.7e\n", Hydro_Con2Eint(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                                                NULL, NULL, CheckMinEint_No, NULL_REAL, Emag_tmp) );
                  }

                  fclose( File );

//                output the failed patch (mainly for recording the sibling information)
#                 ifdef MHD
                  const int MagSg = amr->MagSg[lv];
#                 else
                  const int MagSg = NULL_INT;
#                 endif
#                 ifdef GRAVITY
                  const int PotSg = amr->PotSg[lv];
#                 else
                  const int PotSg = NULL_INT;
#                 endif
                  Output_Patch( lv, PID_Failed, amr->FluSg[lv], MagSg, PotSg, "Unphy" );
               } // if ( ! AutoReduceDt_Continue )
            } // if ( Unphysical(Update, (AutoReduceDt_Continue)?CheckMinEint:CheckMinEtot, Emag_Out) )

            else
            {
//             store the corrected solution
               for (int v=0; v<NCOMP_TOTAL; v++)   h_Flu_Array_F_Out[TID][v][idx_out] = Update[v];


//             replace the original coarse-fine boundary fluxes with the 1st-order fluxes for the flux fix-up
               if ( OPT__FIXUP_FLUX  &&  OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )
               {
                  const int dim_map[3][2] = { {1, 2}, {0, 2}, {0, 1} };

                  int   idx_m, idx_n, idx_store, FaceIdx;
                  bool  Store1stFlux;
                  real *FluxPtr=NULL;

                  for (int d=0; d<3; d++)
                  {
                     Store1stFlux = false;

                     if ( ijk_out[d] == 0  ||  ijk_out[d] == PS1 )
                     {
                        FluxPtr      = FluxL[d];
                        FaceIdx      = 3*d + ijk_out[d]/PS1;
                        Store1stFlux = true;
                     }

                     else if ( ijk_out[d] == PS1-1  ||  ijk_out[d] == PS2-1 )
                     {
                        FluxPtr      = FluxR[d];
                        FaceIdx      = 3*d + (ijk_out[d]+1)/PS1;
                        Store1stFlux = true;
                     }

                     if ( Store1stFlux )
                     {
                        idx_m     = ijk_out[ dim_map[d][1] ];
                        idx_n     = ijk_out[ dim_map[d][0] ];
                        idx_store = idx_m*PS2 + idx_n;

                        for (int v=0; v<NCOMP_TOTAL; v++)   h_Flux_Array[TID][FaceIdx][v][idx_store] = FluxPtr[v];
                     }
                  } // for (int d=0; d<3; d++)
               } // if ( OPT__FIXUP_FLUX  &&  OPT__1ST_FLUX_CORR != FIRST_FLUX_CORR_NONE )

//             record the number of corrected cells
               NCorrThisTime ++;

            } // if ( Unphysical(Update, (AutoReduceDt_Continue)?CheckMinEint:CheckMinEtot, Emag_Out) ) ... else ...
         } // if need correction
      } // i,j,k
   } // for (int TID=0; TID<NPG; TID++)

// "delete" applies to NULL as well
   delete [] Corr1D_InOut;

   } // end of OpenMP parallel region


// operations when CorrectUnphysical() fails
   if ( CorrectUnphy == GAMER_FAILED )
   {
//    if AutoReduceDt_Continue is true, set FluStatus_ThisRank as GAMER_FAILED to rerun Flu_AdvancedDt() with a smaller dt
      if ( AutoReduceDt_Continue )
         FluStatus_ThisRank = GAMER_FAILED;

//    otherwise, terminate the program
      else
         Aux_Error( ERROR_INFO, "fluid solver failed at Rank %d, lv %d, Time %20.14e, Step %ld, Counter %ld ...\n",
                    MPI_Rank, lv, Time[lv], Step, AdvanceCounter[lv] );
   }

// accumulate the total number of corrected cells in one global time-step if CorrectUnphysical() works
   else
   {
      NCorrUnphy[lv] += NCorrThisTime;
   }

} // FUNCTION : CorrectUnphysical



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  StoreElectric
// Description :  Save the coarse-grid electric field on the coarse-fine boundaries for patches at level "lv"
//                --> to be fixed later by CorrectElectric()
//
// Parameter   :  lv          : Target refinement level
//                h_Ele_Array : Host array storing the updated E field
//                NPG         : Number of patch groups to be evaluated
//                PID0_List   : List recording the patch indices with LocalID==0 to be udpated
//                dt          : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void StoreElectric( const int lv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                    const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithElectric )  Aux_Error( ERROR_INFO, "amr->WithElectric is off !!\n" );


#  pragma omp parallel for schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;

//       1. store the E field on faces 0 ~ 5
         for (int s=0; s<6; s++)
         {
//          for bitwise reproducibility, store the E field to be corrected in electric_bitrep[]
            const int EleSize = PS1M1*PS1;
#           ifdef BIT_REP_ELECTRIC
            real (*ElePtr)[EleSize] = ( real (*)[EleSize] )amr->patch[0][lv][PID]->electric_bitrep[s];
#           else
            real (*ElePtr)[EleSize] = ( real (*)[EleSize] )amr->patch[0][lv][PID]->electric[s];
#           endif

            if ( ElePtr != NULL )
            {
               int face_idx, disp_m, disp_n;

               switch ( s )
               {
                  case 0:  case 1:
                     face_idx = TABLE_02( LocalID, 'x', 0, 1 ) + s;
                     disp_m   = TABLE_02( LocalID, 'z', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'y', 0, PS1 );
                     break;

                  case 2:  case 3:
                     face_idx = TABLE_02( LocalID, 'y', 1, 2 ) + s;
                     disp_m   = TABLE_02( LocalID, 'x', 0, PS1 ); // (x,z) instead of (z,x)
                     disp_n   = TABLE_02( LocalID, 'z', 0, PS1 );
                     break;

                  case 4:  case 5:
                     face_idx = TABLE_02( LocalID, 'z', 2, 3 ) + s;
                     disp_m   = TABLE_02( LocalID, 'y', 0, PS1 );
                     disp_n   = TABLE_02( LocalID, 'x', 0, PS1 );
                     break;

                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "s", s );
               }

               const int disp_m_p1 = disp_m + 1;
               const int disp_n_p1 = disp_n + 1;

               for (int m=0; m<PS1M1; m++)   {  const int mm = m + disp_m_p1;
               for (int n=0; n<PS1;   n++)   {  const int nn = n + disp_n;

                  ElePtr[0][ m*PS1   + n ] = dt*h_Ele_Array[TID][face_idx][0][ mm*PS2    + nn ];
               }}

               for (int m=0; m<PS1;   m++)   {  const int mm = m + disp_m;
               for (int n=0; n<PS1M1; n++)   {  const int nn = n + disp_n_p1;

                  ElePtr[1][ m*PS1M1 + n ] = dt*h_Ele_Array[TID][face_idx][1][ mm*PS2P1 + nn ];
               }}
            } // if ( ElePtr != NULL )
         } // for (int s=0; s<6; s++)


//       2. store the E field on edges 6 ~ 17
         real *ElePtr[12];

         for (int s=6; s<18; s++)
         {
//          for bitwise reproducibility, store the E field to be corrected in electric_bitrep[]
#           ifdef BIT_REP_ELECTRIC
            ElePtr[s-6] = amr->patch[0][lv][PID]->electric_bitrep[s];
#           else
            ElePtr[s-6] = amr->patch[0][lv][PID]->electric[s];
#           endif
         }

         if ( ElePtr[0] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x',   0, PS1 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 3, 4 );
            for (int n=0; n<PS1; n++)  ElePtr[ 0][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[1] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 3, 4 );
            for (int n=0; n<PS1; n++)  ElePtr[ 1][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[2] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x',   0, PS1 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 4, 5 );
            for (int n=0; n<PS1; n++)  ElePtr[ 2][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[3] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'x', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'z', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'y', 4, 5 );
            for (int n=0; n<PS1; n++)  ElePtr[ 3][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[4] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y',   0, PS1 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 6, 7 );
            for (int n=0; n<PS1; n++)  ElePtr[ 4][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[5] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 6, 7 );
            for (int n=0; n<PS1; n++)  ElePtr[ 5][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[6] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y',   0, PS1 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 7, 8 );
            for (int n=0; n<PS1; n++)  ElePtr[ 6][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[7] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'y', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'x', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'z', 7, 8 );
            for (int n=0; n<PS1; n++)  ElePtr[ 7][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[8] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z',   0, PS1 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 0, 1 );
            for (int n=0; n<PS1; n++)  ElePtr[ 8][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[9] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 0, 1 );
            for (int n=0; n<PS1; n++)  ElePtr[ 9][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[10] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z',   0, PS1 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 1, 2 );
            for (int n=0; n<PS1; n++)  ElePtr[10][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }
         if ( ElePtr[11] != NULL ) {
            const int disp     = TABLE_02( LocalID, 'z', PS1, PS2 )*PS2 + TABLE_02( LocalID, 'y', 0, PS1 );
            const int face_idx = TABLE_02( LocalID, 'x', 1, 2 );
            for (int n=0; n<PS1; n++)  ElePtr[11][n] = dt*h_Ele_Array[TID][face_idx][0][ n + disp ];
         }

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : StoreElectric



//-------------------------------------------------------------------------------------------------------
// Function    :  CorrectElectric
// Description :  Use the electric field on the coarse-fine boundaries at level SonLv to correct the field
//                at level FaLv
//
// Note        :  1. patch->ele_corrected[] is initialized as false in Flu_AdvanceDt()
//
// Parameter   :  SonLv       : Target refinement level
//                h_Ele_Array : Array storing the updated E field
//                NPG         : Number of patch groups to be evaluated
//                PID0_List   : List recording the patch indices with LocalID==0 to be udpated
//                dt          : Evolution time-step
//-------------------------------------------------------------------------------------------------------
void CorrectElectric( const int SonLv, const real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                      const int NPG, const int *PID0_List, const real dt )
{

// check
   if ( !amr->WithElectric )  Aux_Error( ERROR_INFO, "amr->WithElectric is off !!\n" );


// nothing to do on the root level
   if ( SonLv == 0 )    return;


// 1. correct the E field on faces 0 ~ 5
   const int  FaLv         = SonLv - 1;
   const int  face_idx[6]  = { 0, 2, 3, 5, 6, 8 };
   const int  MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
   const int  didx[2]      = { 1, PS2P1 };
   const int  EleSize      = PS1M1*PS1;
   const real dt_2         = (real)0.5*dt;

#  pragma omp parallel
   {
#     pragma omp for schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         const int PID0  = PID0_List[TID];
         const int FaPID = amr->patch[0][SonLv][PID0]->father;

#        ifdef GAMER_DEBUG
         if ( FaPID < 0 )
            Aux_Error( ERROR_INFO, "FaPID = %d < 0 (SonLv %d, PID0 %d) !!\n", FaPID, SonLv, PID0 );
#        endif

         for (int s=0; s<6; s++)
         {
            const int FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#           ifdef GAMER_DEBUG
            if ( FaSibPID == -1 )
               Aux_Error( ERROR_INFO, "FaSibPID == -1 (FaLv %d, FaPID %d, s %d, PID0 %d) !!\n", FaLv, FaPID, s, PID0 );
#           endif

//          skip patches adjacent to non-periodic boundaries
            if ( FaSibPID < -1 )    continue;

//          for AUTO_REDUCE_DT, store the updated electrid field in the temporary array electric_tmp[] since
//          we may need to abandon them if the fluid solver fails
            real (*CoarseElePtr)[EleSize] = ( real (*)[EleSize] )
                                            (  ( AUTO_REDUCE_DT ) ? amr->patch[0][FaLv][FaSibPID]->electric_tmp[ MirrorSib[s] ] :
                                                                    amr->patch[0][FaLv][FaSibPID]->electric    [ MirrorSib[s] ]  );

//          skip patches not adjacent to coarse-fine boundaries
            if ( CoarseElePtr == NULL )   continue;

//          correct the coarse-grid E field
            for (int t=0; t<EleSize; t++)
            {
                const int n      = t % PS1;
                const int m      = t / PS1 + 1;
                const int idx    = 2*( m*PS2 + n );
                const int idx_p1 = idx + didx[0];

                CoarseElePtr[0][t] -= dt_2*( h_Ele_Array[TID][ face_idx[s] ][0][idx   ] +
                                             h_Ele_Array[TID][ face_idx[s] ][0][idx_p1] );
            }

            for (int t=0; t<EleSize; t++)
            {
                const int n      = t % PS1M1 + 1;
                const int m      = t / PS1M1;
                const int idx    = 2*( m*PS2P1 + n );
                const int idx_p1 = idx + didx[1];

                CoarseElePtr[1][t] -= dt_2*( h_Ele_Array[TID][ face_idx[s] ][1][idx   ] +
                                             h_Ele_Array[TID][ face_idx[s] ][1][idx_p1] );
            }
         } // for (int s=0; s<6; s++)
      } // for (int TID=0; TID<NPG; TID++)

   } // OpenMP parallel region


// 2. correct the E field on edges 6 ~ 17

// group sibling directions into four independent groups to avoid multiple OpenMP threads updating
// the E field at the same patch edge
   const int IndSib[4][3] = { {6, 10, 14}, {7, 11, 15}, {8, 12, 16}, {9, 13, 17} };

   for (int g=0; g<4; g++)
   {
#     pragma omp parallel
      {
         for (int w=0; w<3; w++)
         {
            const int s           = IndSib[g][w];  // 6 ~ 17
            const int e           = s - 6;         // 0 ~ 11 (edge index)
            const int xyz         = e / 4;         // (0,0,0,0,1,1,1,1,2,2,2,2)
            const int face_offset = e % 4 / 2;     // (0,0,1,1,0,0,1,1,0,0,1,1)
            const int m_offset    = e % 2;         // (0,1,0,1,0,1,0,1,0,1,0,1)
            const int face_idx    = ( 3 + xyz*3 )%9 + 2*face_offset;

            int  SibID[3], SibSibID[3];
            real AveFineEle[PS1];

            TABLE_SiblingSharingSameEdge( s, SibID, SibSibID );

#           pragma omp for schedule( runtime )
            for (int TID=0; TID<NPG; TID++)
            {
               const int PID0  = PID0_List[TID];
               const int FaPID = amr->patch[0][SonLv][PID0]->father;

               bool AveFineEleReady = false;

#              ifdef GAMER_DEBUG
               if ( FaPID < 0 )
                  Aux_Error( ERROR_INFO, "FaPID = %d < 0 (SonLv %d, PID0 %d) !!\n", FaPID, SonLv, PID0 );
#              endif

               for (int u=0; u<3; u++)
               {
                  const int FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[ SibID[u] ];

#                 ifdef GAMER_DEBUG
                  if ( FaSibPID == -1 )
                     Aux_Error( ERROR_INFO, "FaSibPID == -1 (FaLv %d, FaPID %d, u %d, PID0 %d) !!\n", FaLv, FaPID, u, PID0 );
#                 endif

//                skip patches adjacent to non-periodic boundaries
                  if ( FaSibPID < -1 )    continue;

//                for AUTO_REDUCE_DT, store the updated electrid field in the temporary array electric_tmp[] since
//                we may need to abandon them if the fluid solver fails
                  real *CoarseElePtr= ( AUTO_REDUCE_DT ) ? amr->patch[0][FaLv][FaSibPID]->electric_tmp[ SibSibID[u] ] :
                                                           amr->patch[0][FaLv][FaSibPID]->electric    [ SibSibID[u] ];

//                skip patches not adjacent to coarse-fine boundaries
                  if ( CoarseElePtr == NULL )   continue;

//                skip the E field arrays that have been corrected already
                  if ( amr->patch[0][FaLv][FaSibPID]->ele_corrected[ SibSibID[u]-6 ] )    continue;

//                prepare the fine-grid E field data if not done yet
                  if ( ! AveFineEleReady )
                  {
                      const int m = m_offset*PS2;

                      for (int t=0; t<PS1; t++)
                      {
                         const int n      = t*2;
                         const int idx    = m*PS2 + n;
                         const int idx_p1 = idx + 1;

                         AveFineEle[t] = dt_2*( h_Ele_Array[TID][face_idx][0][idx   ] +
                                                h_Ele_Array[TID][face_idx][0][idx_p1] );
                      }

                      AveFineEleReady = true;
                  } // if ( ! AveFineEleReady )

//                correct the coarse-grid E field
                  for (int t=0; t<PS1; t++)  CoarseElePtr[t] -= AveFineEle[t];

//                mark that this E field array has been corrected to avoid duplicate correction
//                --> necessary because each coarse patch edge can be overlapped with up to three fine patch edges
                  amr->patch[0][FaLv][FaSibPID]->ele_corrected[ SibSibID[u]-6 ] = true;
               } // for (int u=0; u<3; u++)
            }//for (int TID=0; TID<NPG; TID++)
         } // for (int w=0; w<3; w++)
      } // OpenMP parallel region
   } // for (int g=0; g<4; g++)

} // FUNCTION : CorrectElectric

#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )
