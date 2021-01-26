#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void BC_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                              const int lv, double AuxArray[] );

// this function pointer must be set by a test problem initializer
void (*BC_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                     const int lv, double AuxArray[] ) = NULL;

#ifdef MHD
extern void (*BC_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                   const int lv, double AuxArray[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User_Template
// Description :  User-specified boundary condition template
//
// Note        :  1. Invoked by Flu_BoundaryCondition_User() using the function pointer
//                   "BC_User_Ptr", which must be set by a test problem initializer
//                2. Always return NCOMP_TOTAL variables
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//                   --> Just like the fluid initialization routine
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                       const int lv, double AuxArray[] )
{

// put your B.C. here
// ##########################################################################################################
// Example 1 : set to time-independent values for HYDRO
   /*
   const real Dens0 = 1.0;
   const real Vx    = 1.25e-1;
   const real Vy    = 2.30e-1;
   const real Vz    = 3.70e-1;
   const real Pres0 = 1.0;
   const real Emag0 = 0.0;    // must be zero here even for MHD

   real Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
#  if ( NCOMP_PASSIVE > 0 )
   real Passive[NCOMP_PASSIVE];
#  else
   real *Passive = NULL;
#  endif

   Dens       = Dens0 + 0.2*exp(  -(  SQR(1.1*x-0.5*amr->BoxSize[0])
                                     +SQR(2.2*y-0.5*amr->BoxSize[1])
                                     +SQR(3.3*z-0.5*amr->BoxSize[2]) ) / SQR( 1.8*amr->BoxSize[2] )  );
   MomX       = Dens*Vx;
   MomY       = Dens*Vy;
   MomZ       = Dens*Vz;
   Pres       = Pres0*(  2.0 + sin( 2.0*M_PI*(4.5*x+5.5*y*6.5*z)/amr->BoxSize[2] )  );
#  if ( NCOMP_PASSIVE > 0 )
// Passive[X] = ...;
#  endif
   Eint       = EoS_DensPres2Eint_CPUPtr( Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Etot       = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, Emag0 );

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
#  if ( NCOMP_PASSIVE > 0 )
   fluid[XXXX] = ...;
#  endif
   */


// Example 2 : set to time-dependent values for HYDRO

// ##########################################################################################################

} // FUNCTION : BC_User_Template



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_BoundaryCondition_User
// Description :  Fill up the ghost-zone values by the user-specified boundary condition
//
// Note        :  1. Work for Prepare_PatchData(), InterpolateGhostZone(), Refine(), and LB_Refine_GetNewRealPatchList()
//                2. Function pointers "BC_User_Ptr" and "BC_BField_User_Ptr" must be set by a test problem initializer
//                3. User-defined boundary conditions for the magnetic field are set in MHD_BoundaryCondition_User()
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                NVar_Flu       : Number of fluid variables to be prepared (derived variables are NOT included)
//                ArraySizeX/Y/Z : Size of Array including the ghost zones on each side
//                Idx_Start      : Minimum array indices
//                Idx_End        : Maximum array indices
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                Time           : Current physical time
//                dh             : Cell size
//                Corner         : Physcial coordinates at the center of the cell (0,0,0) --> Array[0]
//                TVar           : Target variables to be prepared --> only used for preparing the derived variables
//                lv             : Refinement level
//
// Return      :  Array
//-------------------------------------------------------------------------------------------------------
void Flu_BoundaryCondition_User( real *Array, const int NVar_Flu, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                                 const int TFluVarIdxList[], const double Time, const double dh, const double *Corner,
                                 const long TVar, const int lv )
{

// check
   if ( BC_User_Ptr == NULL )
      Aux_Error( ERROR_INFO, "BC_User_Ptr == NULL for user-specified boundary conditions !!\n" );

#  ifdef MHD
   if ( BC_BField_User_Ptr == NULL )
      Aux_Error( ERROR_INFO, "BC_BField_User_Ptr == NULL for user-specified boundary conditions !!\n" );
#  endif


   const double x0 = Corner[0] + (double)Idx_Start[0]*dh;   // starting x,y,z coordinates
   const double y0 = Corner[1] + (double)Idx_Start[1]*dh;
   const double z0 = Corner[2] + (double)Idx_Start[2]*dh;

#  if   ( MODEL == HYDRO )
#  ifdef MHD
   const double dh_2             = 0.5*dh;
#  endif
#  ifdef SRHD
   const bool   PrepLrtz         = ( TVar   & _LORENTZ_FACTOR ) ? true : false; // Lorentz factor
#  endif
   const bool   CheckMinPres_Yes = true;
   const bool   PrepVx           = ( TVar & _VELX ) ? true : false;
   const bool   PrepVy           = ( TVar & _VELY ) ? true : false;
   const bool   PrepVz           = ( TVar & _VELZ ) ? true : false;
   const bool   PrepPres         = ( TVar & _PRES ) ? true : false;
   const bool   PrepTemp         = ( TVar & _TEMP ) ? true : false;

#  elif ( MODEL == ELBDM )
// no derived variables yet

#  else
#  error : unsupported MODEL !!
#  endif


// 1D array -> 3D array
   real (*Array3D)[ArraySizeZ][ArraySizeY][ArraySizeX] = ( real (*)[ArraySizeZ][ArraySizeY][ArraySizeX] )Array;


// set the boundary values
   int    i, j, k, v2;
   real   BVal[NCOMP_TOTAL];
   double x, y, z;

   for (k=Idx_Start[2], z=z0; k<=Idx_End[2]; k++, z+=dh)
   for (j=Idx_Start[1], y=y0; j<=Idx_End[1]; j++, y+=dh)
   for (i=Idx_Start[0], x=x0; i<=Idx_End[0]; i++, x+=dh)
   {
//    1. primary variables
//    get the boundary values of all NCOMP_TOTAL fields
      BC_User_Ptr( BVal, x, y, z, Time, lv, NULL );

//    add the magnetic energy for MHD
#     if ( MODEL == HYDRO )
#     ifdef MHD
      real Emag, BxL, BxR, Bx, ByL, ByR, By, BzL, BzR, Bz, B3v[NCOMP_MAG];

      BC_BField_User_Ptr( B3v, x-dh_2, y,      z,      Time, lv, NULL );   BxL = B3v[MAGX];
      BC_BField_User_Ptr( B3v, x+dh_2, y,      z,      Time, lv, NULL );   BxR = B3v[MAGX];
      BC_BField_User_Ptr( B3v, x,      y-dh_2, z,      Time, lv, NULL );   ByL = B3v[MAGY];
      BC_BField_User_Ptr( B3v, x,      y+dh_2, z,      Time, lv, NULL );   ByR = B3v[MAGY];
      BC_BField_User_Ptr( B3v, x,      y,      z-dh_2, Time, lv, NULL );   BzL = B3v[MAGZ];
      BC_BField_User_Ptr( B3v, x,      y,      z+dh_2, Time, lv, NULL );   BzR = B3v[MAGZ];

      Bx   = (real)0.5*( BxL + BxR );
      By   = (real)0.5*( ByL + ByR );
      Bz   = (real)0.5*( BzL + BzR );
      Emag = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );

      BVal[ENGY] += Emag;

#     else
      const real Emag = NULL_REAL;
#     endif
#     endif // #ifdef ( MODEL == HYDRO )

//    store results to the output array
      for (int v=0; v<NVar_Flu; v++)   Array3D[v][k][j][i] = BVal[ TFluVarIdxList[v] ];


//    2. derived variables
      v2 = NVar_Flu;

#     if   ( MODEL == HYDRO )
#     ifdef SRHD
      real LorentzFactor, Prim[NCOMP_FLUID];

      Hydro_Con2Pri( BVal, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                     (real)NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                     EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                     EoS_AuxArray_Int, h_EoS_Table, NULL, &LorentzFactor );

      if ( PrepLrtz )   Array3D[ v2 ++ ][k][j][i] = LorentzFactor;
#     else
      if ( PrepVx   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMX] / BVal[DENS];
      if ( PrepVy   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMY] / BVal[DENS];
      if ( PrepVz   )   Array3D[ v2 ++ ][k][j][i] = BVal[MOMZ] / BVal[DENS];
#     endif
      if ( PrepPres )   Array3D[ v2 ++ ][k][j][i] = Hydro_Con2Pres( BVal[DENS], BVal[MOMX], BVal[MOMY],
                                                                    BVal[MOMZ], BVal[ENGY], BVal+NCOMP_FLUID,
                                                                    CheckMinPres_Yes, MIN_PRES, Emag,
                                                                    EoS_DensEint2Pres_CPUPtr,
                                                                    EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                                    h_EoS_Table, NULL );
      if ( PrepTemp )   Array3D[ v2 ++ ][k][j][i] = Hydro_Con2Temp( BVal[DENS], BVal[MOMX], BVal[MOMY],
                                                                    BVal[MOMZ], BVal[ENGY], BVal+NCOMP_FLUID,
                                                                    CheckMinPres_Yes, MIN_PRES, Emag,
                                                                    EoS_DensEint2Pres_CPUPtr,
                                                                    EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

#     elif ( MODEL == ELBDM )
//    no derived variables yet

#     else
#     error : unsupported MODEL !!
#     endif
   } // k,j,i

} // FUNCTION : Flu_BoundaryCondition_User
