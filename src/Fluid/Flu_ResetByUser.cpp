#include "GAMER.h"

// declare as static so that other functions cannot invoke them directly and must use the function pointers
static bool Flu_ResetByUser_Func_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                           const int lv, double AuxArray[] );
static void Flu_ResetByUser_API_Default( const int lv, const int FluSg, const double TTime );

// this function pointer **must** be set by a test problem initializer
bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] ) = NULL;
// this function pointer **may** be overwritten by a test problem initializer
void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime ) = Flu_ResetByUser_API_Default;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func_Template
// Description :  Function template to reset the fluid field
//
// Note        :  1. Invoked by Flu_ResetByUser_API_Default() and Model_Init_ByFunction_AssignData() using the
//                   function pointer "Flu_ResetByUser_Func_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked when constructing the initial condition in
//                   Model_Init_ByFunction_AssignData() and after each update in Flu_ResetByUser_API_Default()
//                3. Input fluid[] stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in Flu_ResetByUser_API_Default() and
//                       Model_Init_ByFunction_AssignData()
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (original) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_Func_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                    const int lv, double AuxArray[] )
{

// Example : reset fluid variables to extremely small values if the cell is within a specific sphere
   /*
   const real dr[3]   = { x-amr->BoxCenter[0], y-amr->BoxCenter[1], z-amr->BoxCenter[2] };
   const real r       = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   const real TRad                      = 0.3;
   const real MinDens                   = 1.0e-10;
#  if ( NCOMP_PASSIVE > 0 )
   const real MinPassive[NCOMP_PASSIVE] = { ... };
#  else
   const real *MinPassive = NULL;
#  endif
   const real MinPres                   = 1.0e-10;
   const real MinEint                   = EoS_DensPres2Eint_CPUPtr( MinDens, MinPres, MinPassive,
                                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   const real MinEmag                   = 0.0;  // assuming MHD is not adopted

   if ( r <= TRad )
   {
//    set active scalars
      fluid[DENS] = MinDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], MinEint, MinEmag );

//    set passive scalars
#     if ( NCOMP_PASSIVE > 0 )
//    fluid[XXXX] = ...;
#     endif

      return true;
   }

   else
      return false;
   */

   return false;

} // FUNCTION : Flu_ResetByUser_Func_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API_Default
// Description :  Default API for resetting the fluid array
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by Flu_AdvanceDt(), Gra_AdvanceDt(), and Grackle_AdvanceDt() using the
//                   function pointer "Flu_ResetByUser_API_Ptr"
//                   --> This function pointer may be reset by a test problem initializer, in which case
//                       this funtion will become useless
//                3. Currently NOT applied to the input uniform array
//                   --> Init_ByFile() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                TTime : Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API_Default( const int lv, const int FluSg, const double TTime )
{

// check
   if ( Flu_ResetByUser_Func_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL for OPT__RESET_FLUID !!\n" );


   const double dh = amr->dh[lv];

   bool   Reset;
   real   fluid[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


#  pragma omp parallel for private( Reset, fluid, x, y, z, x0, y0, z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

      for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;
      for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;
      for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//       reset this cell
         Reset = Flu_ResetByUser_Func_Ptr( fluid, x, y, z, TTime, lv, NULL );

//       operations necessary only when this cell has been reset
//       --> SRHD do not check and remedy unphysically reseted cells temporarily.
         if ( Reset )
         {
#           if ( MODEL == HYDRO && !defined SRHD )
#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#           else
            const real Emag = NULL_REAL;
#           endif

//          apply density and internal energy floors
            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
            fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                    MIN_EINT, Emag );

//          calculate the dual-energy variable (entropy or internal energy)
#           if   ( DUAL_ENERGY == DE_ENPY )
            fluid[ENPY] = Hydro_Con2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Emag,
                                             EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif
#           endif // if ( MODEL == HYDRO )

//          store the reset values
            for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Flu_ResetByUser_API_Default
