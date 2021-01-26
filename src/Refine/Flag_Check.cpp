#include "GAMER.h"

static bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold );
static bool Check_Curl( const int i, const int j, const int k,
                        const real vx[][PS1][PS1], const real vy[][PS1][PS1], const real vz[][PS1][PS1],
                        const double Threshold );
extern bool (*Flag_Region_Ptr)( const int i, const int j, const int k, const int lv, const int PID );
extern bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Check
// Description :  Check if the target cell (i,j,k) satisfies the refinement criteria
//
// Note        :  1. Useless input arrays are set to NULL (e.g, Pot[] if GRAVITY is off)
//                2. For OPT__FLAG_USER, the function pointer "Flag_User_Ptr" must be set by a
//                   test problem initializer
//
// Parameter   :  lv           : Target refinement level
//                PID          : Target patch ID
//                i,j,k        : Indices of the target cell
//                dv           : Cell volume at the target level
//                Fluid        : Input fluid array (with NCOMP_TOTAL components)
//                Pot          : Input potential array
//                MagCC        : Input cell-centered B field array
//                Vel          : Input velocity array
//                Pres         : Input pressure array
//                Lohner_Ave   : Input array storing the averages for the Lohner error estimator
//                Lohner_Slope : Input array storing the slopes for the Lohner error estimator
//                Lohner_NVar  : Number of variables stored in Lohner_Ave and Lohner_Slope
//                ParCount     : Input array storing the number of particles on each cell
//                               (note that it has the **real** type)
//                ParDens      : Input array storing the particle mass density on each cell
//                JeansCoeff   : Pi*GAMMA/(SafetyFactor^2*G), where SafetyFactor = FlagTable_Jeans[lv]
//                               --> Flag if dh^2 > JeansCoeff*Pres/Dens^2
//
// Return      :  "true"  if any  of the refinement criteria is satisfied
//                "false" if none of the refinement criteria is satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, const real dv,
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                 const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1], const real LorentzFactor[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                 const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1], const real JeansCoeff )
{

   bool Flag = false;

// check whether the input cell is within the regions allowed to be refined
// ===========================================================================================
   if ( OPT__FLAG_REGION )
   {
      if (  !Flag_Region_Ptr( i, j, k, lv, PID )  )    return false;
   }


#  ifdef PARTICLE
// check the number of particles on each cell
// ===========================================================================================
   if ( OPT__FLAG_NPAR_CELL )
   {
      Flag |= ( ParCount[k][j][i] > FlagTable_NParCell[lv] );
      if ( Flag )    return Flag;
   }


// check the particle mass on each cell
// ===========================================================================================
   if ( OPT__FLAG_PAR_MASS_CELL )
   {
      Flag |= ( ParDens[k][j][i]*dv > FlagTable_ParMassCell[lv] );
      if ( Flag )    return Flag;
   }
#  endif


#  ifdef DENS
// check density magnitude
// ===========================================================================================
   if ( OPT__FLAG_RHO )
   {
      Flag |= ( Fluid[DENS][k][j][i] > FlagTable_Rho[lv] );
      if ( Flag )    return Flag;
   }


// check density gradient
// ===========================================================================================
   if ( OPT__FLAG_RHO_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Fluid[DENS][0][0][0], FlagTable_RhoGradient[lv] );
      if ( Flag )    return Flag;
   }
#  endif


#  if ( MODEL == HYDRO )
// check pressure gradient
// ===========================================================================================
   if ( OPT__FLAG_PRES_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Pres[0][0][0], FlagTable_PresGradient[lv] );
      if ( Flag )    return Flag;
   }


// check gradient of reduced energy density
// ===========================================================================================
   if ( OPT__FLAG_ENGY_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &Fluid[ENGY][0][0][0], FlagTable_EngyGradient[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check vorticity
// ===========================================================================================
#  if ( MODEL == HYDRO  &&  !defined SRHD )
   if ( OPT__FLAG_VORTICITY )
   {
      Flag |= Check_Curl( i, j, k, Vel[0], Vel[1], Vel[2], FlagTable_Vorticity[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check current density in MHD
// ===========================================================================================
#  ifdef MHD
   if ( OPT__FLAG_CURRENT )
   {
      Flag |= Check_Curl( i, j, k, MagCC[0], MagCC[1], MagCC[2], FlagTable_Current[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check M/D in SRHD
// ===========================================================================================
#  ifdef SRHD
   if ( OPT__FLAG_MOM_OVER_DENS )
   {
      real Cons[NCOMP_FLUID], M_D;
      Cons[DENS]=Fluid[DENS][k][j][i];
      Cons[MOMX]=Fluid[MOMX][k][j][i];
      Cons[MOMY]=Fluid[MOMY][k][j][i];
      Cons[MOMZ]=Fluid[MOMZ][k][j][i];

      M_D = FABS( Cons[MOMX] ) + FABS( Cons[MOMY] ) + FABS( Cons[MOMZ] );

      M_D /= Cons[DENS];

      Flag |= ( M_D > FlagTable_Mom_Over_Dens[lv] );
      if ( Flag )    return Flag;
   }


// check magnitude of 4-velocity
// ===========================================================================================
   if ( OPT__FLAG_4VELOCITY )
   {
      real Dens, MomX, MomY, MomZ, Engy;
      real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID], U, LorentzFactor;

      Cons[DENS]=Fluid[DENS][k][j][i];
      Cons[MOMX]=Fluid[MOMX][k][j][i];
      Cons[MOMY]=Fluid[MOMY][k][j][i];
      Cons[MOMZ]=Fluid[MOMZ][k][j][i];
      Cons[ENGY]=Fluid[ENGY][k][j][i];

      Hydro_Con2Pri( Cons, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                     (real)NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                     EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                     EoS_AuxArray_Int, h_EoS_Table, NULL, &LorentzFactor );

      U = SQRT( LorentzFactor*LorentzFactor - (real)1.0 );

      Flag |= ( U > FlagTable_4Velocity[lv] );

      if ( Flag )    return Flag;
   }


// check gradient of Lorentz factor
// ===========================================================================================
   if ( OPT__FLAG_LORENTZ_GRADIENT )
   {
      Flag |= Check_Gradient( i, j, k, &LorentzFactor[0][0][0], FlagTable_LorentzFactorGradient[lv] );
      if ( Flag )    return Flag;
   }
#  endif


// check Jeans length
// ===========================================================================================
#  if ( MODEL == HYDRO  &&  defined GRAVITY )
   if ( OPT__FLAG_JEANS )
   {
#     ifdef GAMER_DEBUG
      if ( Pres == NULL )  Aux_Error( ERROR_INFO, "Pres == NULL !!\n" );
#     endif

      const real Dens_1Cell = Fluid[DENS][k][j][i];
      const real Pres_1Cell = Pres       [k][j][i];

      Flag |= (  SQR(amr->dh[lv]) > JeansCoeff*Pres_1Cell/SQR( Dens_1Cell )  );
      if ( Flag )    return Flag;
   }
#  endif


// check ELBDM energy density
// ===========================================================================================
#  if ( MODEL == ELBDM )
   if ( OPT__FLAG_ENGY_DENSITY )
   {
      Flag |= ELBDM_Flag_EngyDensity( i, j, k, &Fluid[REAL][0][0][0], &Fluid[IMAG][0][0][0],
                                      FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
      if ( Flag )    return Flag;
   }
#  endif


// check Lohner's error estimator
// ===========================================================================================
   if ( Lohner_NVar > 0 )
   {
//    check Lohner only if density is greater than the minimum threshold
#     ifdef DENS
      if ( Fluid[DENS][k][j][i] >= FlagTable_Lohner[lv][3] )
#     endif
      Flag |= Flag_Lohner( i, j, k, OPT__FLAG_LOHNER_FORM, Lohner_Var, Lohner_Ave, Lohner_Slope, Lohner_NVar,
                           FlagTable_Lohner[lv][0], FlagTable_Lohner[lv][1], FlagTable_Lohner[lv][2] );
      if ( Flag )    return Flag;
   }


// check user-defined criteria
// ===========================================================================================
   if ( OPT__FLAG_USER )
   {
      if ( Flag_User_Ptr != NULL )
      {
         Flag |= Flag_User_Ptr( i, j, k, lv, PID, FlagTable_User[lv] );
         if ( Flag )    return Flag;
      }

      else
         Aux_Error( ERROR_INFO, "Flag_User_Ptr == NULL for OPT__FLAG_USER !!\n" );
   }


   return Flag;

} // FUNCTION : Flag_Check



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Gradient
// Description :  Check if the gradient of the input data at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Size of the array "Input" should be PATCH_SIZE^3
//                2. For cells adjacent to the patch boundaries, only first-order approximation is adopted
//                   to estimate gradient. Otherwise, second-order approximation is adopted.
//                   --> Advantage: NO need to prepare the ghost-zone data for the target patch
//
// Parameter   :  i,j,k     : Indices of the target cell in the array "Input"
//                Input     : Input array
//                Threshold : Threshold for the flag operation
//
// Return      :  "true"  if the gradient is larger           than the given threshold
//                "false" if the gradient is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool Check_Gradient( const int i, const int j, const int k, const real Input[], const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif


   const int ijk[3]  = { i, j, k };
   const int Idx     = k*PS1*PS1 + j*PS1 + i;
   const int dIdx[3] = { 1, PS1, PS1*PS1 };

   int  Idx_p, Idx_m;
   real _dh, Self, Gradient;
   bool Flag = false;

   for (int d=0; d<3; d++)
   {
      switch ( ijk[d] )
      {
         case 0     : Idx_m = Idx;           Idx_p = Idx+dIdx[d];    _dh = (real)1.0;  break;
         case PS1-1 : Idx_m = Idx-dIdx[d];   Idx_p = Idx;            _dh = (real)1.0;  break;
         default    : Idx_m = Idx-dIdx[d];   Idx_p = Idx+dIdx[d];    _dh = (real)0.5;  break;
      }

      Self     = Input[Idx];
      Gradient = _dh*( Input[Idx_p] - Input[Idx_m] );
      Flag    |= (  FABS( Gradient/Self ) > Threshold  );

      if ( Flag )    return Flag;
   } // for (int d=0; d<3; d++)

   return Flag;

} // FUNCTION : Check_Gradient



//-------------------------------------------------------------------------------------------------------
// Function    :  Check_Curl
// Description :  Check if the curl of the input vector at the cell (i,j,k) exceeds the given threshold
//
// Note        :  1. Flag if |curl(v)|*dh/|v| > threshold
//                2. For cells adjacent to the patch boundaries, only first-order approximation is adopted
//                   to estimate derivatives. Otherwise, second-order approximation is adopted.
//                   --> Advantage: NO need to prepare the ghost-zone data for the target patch
//                3. Size of the input arrays "vx/y/z" should be PATCH_SIZE^3
//                   --> They should store **cell-centered** values
//
// Parameter   :  i,j,k     : Target array indices
//                vx/y/z    : Input vectors
//                Threshold : Refinement threshold
//
// Return      :  "true"  if |curl(v)|*dh/|v| >  threshold
//                "false" if |curl(v)|*dh/|v| <= threshold
//-------------------------------------------------------------------------------------------------------
bool Check_Curl( const int i, const int j, const int k,
                 const real vx[][PS1][PS1], const real vy[][PS1][PS1], const real vz[][PS1][PS1],
                 const double Threshold )
{

// check
#  ifdef GAMER_DEBUG
   if (  i < 0  ||  i >= PS1  ||  j < 0 ||  j >= PS1  ||  k < 0  ||  k >= PS1   )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );
#  endif

   int  im, ip, jm, jp, km, kp;
   real _dx, _dy, _dz, v2, wx, wy, wz, w2;
   bool Flag = false;


// check if the target cell indices are adjacent to the patch boundaries
   if      ( i == 0     ) { im = 0;       ip = 1;       _dx = (real)1.0; }
   else if ( i == PS1-1 ) { im = PS1-2;   ip = PS1-1;   _dx = (real)1.0; }
   else                   { im = i-1;     ip = i+1;     _dx = (real)0.5; }

   if      ( j == 0     ) { jm = 0;       jp = 1;       _dy = (real)1.0; }
   else if ( j == PS1-1 ) { jm = PS1-2;   jp = PS1-1;   _dy = (real)1.0; }
   else                   { jm = j-1;     jp = j+1;     _dy = (real)0.5; }

   if      ( k == 0     ) { km = 0;       kp = 1;       _dz = (real)1.0; }
   else if ( k == PS1-1 ) { km = PS1-2;   kp = PS1-1;   _dz = (real)1.0; }
   else                   { km = k-1;     kp = k+1;     _dz = (real)0.5; }


// calculate magnitude
   v2 = SQR( vx[k][j][i] ) + SQR( vy[k][j][i] ) + SQR( vz[k][j][i] );


// calculate w=curl(v)*dh
   wx = _dy*( vz[k ][jp][i ] - vz[k ][jm][i ] ) - _dz*( vy[kp][j ][i ] - vy[km][j ][i ] );
   wy = _dz*( vx[kp][j ][i ] - vx[km][j ][i ] ) - _dx*( vz[k ][j ][ip] - vz[k ][j ][im] );
   wz = _dx*( vy[k ][j ][ip] - vy[k ][j ][im] ) - _dy*( vx[k ][jp][i ] - vx[k ][jm][i ] );
   w2 = SQR(wx) + SQR(wy) + SQR(wz);


// flag if |curl(v)|*dh/|v| > threshold
   Flag = ( w2/v2 > SQR(Threshold) );

   return Flag;

} // FUNCTION : Check_Curl
