#ifndef __CUFLU_RIEMANNSOLVER_HLLE__
#define __CUFLU_RIEMANNSOLVER_HLLE__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real* const PresIn );

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLE
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : (a) Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                         (c) Einfeldt et al., J. Comput. Phys., 92, 273 (1991)
//                3. Wave-speed estimator is set by HLLE_WAVESPEED in CUFLU.h
//                4. Support general EoS
//                5. Shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out          : Array to store the output flux
//                L/R_In            : Input left/right states (conserved variables)
//                MinDens/Pres      : Density and pressure floors
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] )
{

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL_PLUS_MAG], R[NCOMP_TOTAL_PLUS_MAG];

   for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &L[0], NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &R[0], NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
#  endif


// 2. estimate the maximum wave speeds
// 2-1. compute the left/right states
   const real ZERO             = (real)0.0;
   const real ONE              = (real)1.0;
   const real _TWO             = (real)0.5;
   const bool CheckMinPres_Yes = true;

   real _RhoL, _RhoR, u_L, u_R, Emag_L, Emag_R, P_L, P_R, a2_L, a2_R, Cf_L, Cf_R, MaxV_L, MaxV_R;
#  ifdef MHD
   real Bx_L, By_L, Bz_L, Bx_R, By_R, Bz_R, Bx2_L, Bt2_L, Bx2_R, Bt2_R, B2_L, B2_R;
   real Cax2_L, Cat2_L, Ca2_plus_a2_L, Ca2_min_a2_L, Cf2_min_Cs2_L, Cf2_L;    // "plus"="+", "min"="-", Cs=slow wave
   real Cax2_R, Cat2_R, Ca2_plus_a2_R, Ca2_min_a2_R, Cf2_min_Cs2_R, Cf2_R;
#  endif

   _RhoL = ONE / L[0];
   _RhoR = ONE / R[0];
   u_L   = _RhoL*L[1];
   u_R   = _RhoR*R[1];

#  ifdef MHD
   Bx_L   = L[ MAG_OFFSET + 0 ];
   By_L   = L[ MAG_OFFSET + 1 ];
   Bz_L   = L[ MAG_OFFSET + 2 ];
   Bx_R   = R[ MAG_OFFSET + 0 ];
   By_R   = R[ MAG_OFFSET + 1 ];
   Bz_R   = R[ MAG_OFFSET + 2 ];
   Bx2_L  = SQR( Bx_L );
   Bx2_R  = SQR( Bx_R );
   Bt2_L  = SQR( By_L ) + SQR( Bz_L );
   Bt2_R  = SQR( By_R ) + SQR( Bz_R );
   B2_L   = Bx2_L + Bt2_L;
   B2_R   = Bx2_R + Bt2_R;
   Emag_L = _TWO*B2_L;
   Emag_R = _TWO*B2_R;
#  else
   Emag_L = NULL_REAL;
   Emag_R = NULL_REAL;
#  endif

   P_L   = Hydro_Con2Pres( L[0], L[1], L[2], L[3], L[4], L+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag_L,
                           EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   P_R   = Hydro_Con2Pres( R[0], R[1], R[2], R[3], R[4], R+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag_R,
                           EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   a2_L  = EoS_DensPres2CSqr( L[0], P_L, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   a2_R  = EoS_DensPres2CSqr( R[0], P_R, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &P_L, NULL, "pressure", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &P_R, NULL, "pressure", __FILE__, __FUNCTION__, __LINE__, true );
#  endif

#  if ( defined GAMER_DEBUG  &&  defined MHD )
// longitudinal B field in the left and right states should be the same
   if ( Bx_L != Bx_R )
      printf( "ERROR : Bx_L (%24.17e) != Bx_R (%24.17e) for XYZ %d at file <%s>, line <%d>, function <%s>!!\n",
              Bx_L, Bx_R, XYZ, __FILE__, __LINE__, __FUNCTION__ );
#  endif

// fast wave speed (Cf)
#  ifdef MHD
// left state
   Cax2_L        = Bx2_L*_RhoL;
   Cat2_L        = Bt2_L*_RhoL;
   Ca2_plus_a2_L = Cat2_L + Cax2_L + a2_L;
   Ca2_min_a2_L  = Cat2_L + Cax2_L - a2_L;
   Cf2_min_Cs2_L = SQRT( SQR(Ca2_min_a2_L) + (real)4.0*a2_L*Cat2_L );

   if ( Cat2_L == ZERO )
   {
      if ( Cax2_L >= a2_L )   Cf2_L = Cax2_L;
      else                    Cf2_L = a2_L;
   }

   else
   {
      if ( Cax2_L == ZERO )   Cf2_L = a2_L + Cat2_L;
      else                    Cf2_L = _TWO*( Ca2_plus_a2_L + Cf2_min_Cs2_L );
   } // if ( Cat2_L == ZERO ) ... else ...

   Cf_L = SQRT( Cf2_L );   // Cf2_L is positive definite using the above formula

// right state
   Cax2_R        = Bx2_R*_RhoR;
   Cat2_R        = Bt2_R*_RhoR;
   Ca2_plus_a2_R = Cat2_R + Cax2_R + a2_R;
   Ca2_min_a2_R  = Cat2_R + Cax2_R - a2_R;
   Cf2_min_Cs2_R = SQRT( SQR(Ca2_min_a2_R) + (real)4.0*a2_R*Cat2_R );

   if ( Cat2_R == ZERO )
   {
      if ( Cax2_R >= a2_R )   Cf2_R = Cax2_R;
      else                    Cf2_R = a2_R;
   }

   else
   {
      if ( Cax2_R == ZERO )   Cf2_R = a2_R + Cat2_R;
      else                    Cf2_R = _TWO*( Ca2_plus_a2_R + Cf2_min_Cs2_R );
   } // if ( Cat2_R == ZERO ) ... else ...

   Cf_R = SQRT( Cf2_R );   // Cf2_R is positive definite using the above formula

#  else // #ifdef MHD

   Cf_L = SQRT( a2_L );
   Cf_R = SQRT( a2_R );

#  endif // #ifdef MHD ... else ...


// 2-2a. use the Roe average eigenvalues
#  if   ( HLLE_WAVESPEED == HLL_WAVESPEED_ROE )

// Roe averages
   real H_L, H_R, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real Rho_Roe, _Rho_Roe, u_Roe, v_Roe, w_Roe, H_Roe;
#  ifdef MHD
   real Bx_Roe, By_Roe, Bz_Roe, Bt2_Roe, B2_Roe, B2_Rho_Roe;
#  endif

   H_L   = ( L[4] + P_L )*_RhoL;
   H_R   = ( R[4] + P_R )*_RhoR;
#  ifdef MHD
   H_L  += _TWO*B2_L*_RhoL;
   H_R  += _TWO*B2_R*_RhoR;
#  endif

   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   _RhoL_sqrt      = ONE/RhoL_sqrt;
   _RhoR_sqrt      = ONE/RhoR_sqrt;
   _RhoLR_sqrt_sum = ONE/(RhoL_sqrt + RhoR_sqrt);

   Rho_Roe  = RhoL_sqrt*RhoR_sqrt;
   _Rho_Roe = ONE/Rho_Roe;
   u_Roe    = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v_Roe    = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w_Roe    = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   H_Roe    = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

#  ifdef MHD
   Bx_Roe     = Bx_L;
   By_Roe     = _RhoLR_sqrt_sum*( RhoL_sqrt*By_R + RhoR_sqrt*By_L );
   Bz_Roe     = _RhoLR_sqrt_sum*( RhoL_sqrt*Bz_R + RhoR_sqrt*Bz_L );
   Bt2_Roe    = SQR( By_Roe ) + SQR( Bz_Roe );
   B2_Roe     = SQR( Bx_Roe ) + Bt2_Roe;
   B2_Rho_Roe = B2_Roe*_Rho_Roe;
#  endif

// fast wave speed (Cf_Roe)
//###NOTE: we have assumed a constant-gamma EoS here
// --> otherwise, one needs to specify how to convert (H-0.5*V2, Rho) to Cs^2
// --> see Eq. [A4] in Coleman 2020
#  if ( EOS != EOS_GAMMA )
#     error : ERROR : HLL_WAVESPEED_ROE only works with EOS_GAMMA !!
#  endif

   const real  Gamma    = (real)EoS_AuxArray_Flt[0];
   const real  Gamma_m1 = (real)EoS_AuxArray_Flt[1];
   const real _Gamma    = (real)EoS_AuxArray_Flt[3];
#  ifdef MHD
   const real  Gamma_m2 = Gamma - (real)2.0;
#  endif

   real V2_Roe, a2_Roe, Cf2_Roe, Cf_Roe;
#  ifdef MHD
   real X, Y;  // Eqs. (B15) and (B16) in ref-b
   real Cax2_Roe, Cat2_Roe, Ca2_plus_a2_Roe, Ca2_min_a2_Roe, Cf2_min_Cs2_Roe;
#  endif

   V2_Roe  = SQR( u_Roe ) + SQR( v_Roe ) + SQR( w_Roe );
   a2_Roe  = Gamma_m1*( H_Roe - _TWO*V2_Roe );
#  ifdef MHD
   a2_Roe -= Gamma_m1*B2_Rho_Roe;   // H = 0.5*v^2 + B^2/rho + gamma/(gamma-1)*P/rho
#  endif
   a2_Roe  = Gamma*_Rho_Roe*Hydro_CheckMinPres( a2_Roe*Rho_Roe*_Gamma, MinPres );   // apply pressure floor

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &a2_Roe, NULL, "a2_Roe", __FILE__, __FUNCTION__, __LINE__, true );
#  endif // #ifdef CHECK_UNPHYSICAL_IN_FLUID

#  ifdef MHD
   X               = _TWO*( SQR(By_R-By_L) + SQR(Bz_R-Bz_L) )*SQR( _RhoLR_sqrt_sum );
   X              *= Gamma_m2;
#  ifdef EULERY
   Y               = _TWO*( L[0] + R[0] )*_Rho_Roe;
#  else
   Y               = ONE;
#  endif // #ifdef EULER ... else ...
   Y              *= Gamma_m2;
   a2_Roe         -= X;
   Cax2_Roe        = SQR(Bx_Roe)*_Rho_Roe;
   Cat2_Roe        = ( Gamma_m1 - Y )*Bt2_Roe*_Rho_Roe;
   Ca2_plus_a2_Roe = Cat2_Roe + Cax2_Roe + a2_Roe;
   Ca2_min_a2_Roe  = Cat2_Roe + Cax2_Roe - a2_Roe;
   Cf2_min_Cs2_Roe = SQRT( SQR(Ca2_min_a2_Roe) + (real)4.0*a2_Roe*Cat2_Roe );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &a2_Roe, NULL, "a2_Roe", __FILE__, __FUNCTION__, __LINE__, true );

#  ifdef MHD
   Hydro_CheckUnphysical( NULL, NULL, &Gamma_m1-Y, NULL, "Gamma_m1-Y", __FILE__, __FUNCTION__, __LINE__, true );
#  endif
#  endif // #ifdef CHECK_UNPHYSICAL_IN_FLUID

   if ( Cat2_Roe == ZERO )
   {
      if      ( Cax2_Roe == a2_Roe )   Cf2_Roe = a2_Roe;
      else if ( Cax2_Roe >  a2_Roe )   Cf2_Roe = Cax2_Roe;
      else                             Cf2_Roe = a2_Roe;
   }

   else
   {
      if      ( Cax2_Roe == ZERO )     Cf2_Roe = a2_Roe + Cat2_Roe;
      else                             Cf2_Roe = _TWO*( Ca2_plus_a2_Roe + Cf2_min_Cs2_Roe );
   } // if ( Cat2_Roe == ZERO ) ... else ...

#  else  // #ifdef MHD

   Cf2_Roe = a2_Roe;

#  endif // #ifdef MHD ... else ...

   Cf_Roe = SQRT( Cf2_Roe );

// maximum and minimum eigenvalues
   const real EVal_min = u_Roe - Cf_Roe;
   const real EVal_max = u_Roe + Cf_Roe;

// left/right maximum wave speeds
   MaxV_L = FMIN( EVal_min, u_L-Cf_L );
   MaxV_R = FMAX( EVal_max, u_R+Cf_R );
   MaxV_L = FMIN( MaxV_L, ZERO );
   MaxV_R = FMAX( MaxV_R, ZERO );


// 2-2b. use the primitive variable Riemann solver (PVRS)
#  elif ( HLLE_WAVESPEED == HLL_WAVESPEED_PVRS )

#  ifdef MHD
#     error : HLL_WAVESPEED_PVRS does not support MHD !!
#  endif

// As=a=sound speed in PVRS
   real Rho_PVRS, As_PVRS, RhoAs_PVRS, P_PVRS, Gamma_SL, Gamma_SR, q_L, q_R;

   Rho_PVRS    = _TWO*( L[0] + R[0] );
   As_PVRS     = _TWO*( Cf_L + Cf_R );    // Cf=As for hydro
   RhoAs_PVRS  = Rho_PVRS * As_PVRS;
   P_PVRS      = _TWO*(  ( P_L + P_R ) + ( u_L - u_R )*RhoAs_PVRS  );
   P_PVRS      = Hydro_CheckMinPres( P_PVRS, MinPres );

// for EOS_GAMMA/EOS_ISOTHERMAL, the calculations of Gamma_SL/R can be greatly simplified
// --> results should be exactly the same except for round-off errors
#  if   ( EOS == EOS_GAMMA )
   Gamma_SL    = (real)EoS_AuxArray_Flt[0];
   Gamma_SR    = (real)EoS_AuxArray_Flt[0];
#  elif ( EOS == EOS_ISOTHERMAL )
   Gamma_SL    = ONE;
   Gamma_SR    = ONE;
#  else
   real u_PVRS, Rho_As_PVRS, Rho_SL, Rho_SR, _P;

   u_PVRS      = _TWO*(  ( u_L + u_R ) + ( P_L - P_R )/RhoAs_PVRS  );
   Rho_As_PVRS = Rho_PVRS / As_PVRS;
   Rho_SL      = L[0] + ( u_L - u_PVRS )*Rho_As_PVRS;
   Rho_SR      = R[0] + ( u_PVRS - u_R )*Rho_As_PVRS;
   Rho_SL      = FMAX( Rho_SL, MinDens );
   Rho_SR      = FMAX( Rho_SR, MinDens );
   _P          = ONE / P_PVRS;
// see Eq. [9.8] in Toro 1999 for passive scalars
   Gamma_SL    = EoS_DensPres2CSqr( Rho_SL, P_PVRS, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL )*Rho_SL*_P;
   Gamma_SR    = EoS_DensPres2CSqr( Rho_SR, P_PVRS, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL )*Rho_SR*_P;
#  endif // EOS

   q_L    = ( P_PVRS <= P_L ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SL + ONE )/Gamma_SL*( P_PVRS/P_L - ONE )  );
   q_R    = ( P_PVRS <= P_R ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SR + ONE )/Gamma_SR*( P_PVRS/P_R - ONE )  );
   MaxV_L = u_L - Cf_L*q_L;   // Cf=As for hydro
   MaxV_R = u_R + Cf_R*q_R;
   MaxV_L = FMIN( MaxV_L, ZERO );
   MaxV_R = FMAX( MaxV_R, ZERO );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &q_L, NULL, "q_L", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &q_R, NULL, "q_R", __FILE__, __FUNCTION__, __LINE__, true );
#  endif


// 2-2c. use the min/max of the left and right eigenvalues
#  elif ( HLLE_WAVESPEED == HLL_WAVESPEED_DAVIS )
   MaxV_L = FMIN( u_L-Cf_L, u_R-Cf_R );
   MaxV_R = FMAX( u_L+Cf_L, u_R+Cf_R );
   MaxV_L = FMIN( MaxV_L, ZERO );
   MaxV_R = FMAX( MaxV_R, ZERO );


#  else
#  error : ERROR : unsupported HLLE_WAVESPEED !!
#  endif // HLLE_WAVESPEED


// 3. evaluate the left and right fluxes along the maximum wave speeds
#  ifdef MHD
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4, MAG_OFFSET+1, MAG_OFFSET+2 };
#  else
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };
#  endif
   real Flux_L[NCOMP_TOTAL_PLUS_MAG], Flux_R[NCOMP_TOTAL_PLUS_MAG];  // use NCOMP_TOTAL_PLUS_MAG for Hydro_Con2Flux()

   Hydro_Con2Flux( 0, Flux_L, L, MinPres, NULL, NULL, NULL, NULL, &P_L );
   Hydro_Con2Flux( 0, Flux_R, R, MinPres, NULL, NULL, NULL, NULL, &P_R );

   for (int v=0; v<NWAVE; v++)
   {
      Flux_L[ idx_wave[v] ] -= MaxV_L*L[ idx_wave[v] ];
      Flux_R[ idx_wave[v] ] -= MaxV_R*R[ idx_wave[v] ];
   }


// 4. evaluate the HLLE fluxes
// deal with the special case of MaxV_L=MaxV_R=0
//###REVISE: should it return zero flux due to symmetry?
   if ( MaxV_L == ZERO  &&  MaxV_R == ZERO )
   {
      for (int v=0; v<NWAVE; v++)
         Flux_Out[ idx_wave[v] ] = Flux_L[ idx_wave[v] ];   // assuming Flux_L=Flux_R
   }

   else
   {
      const real _MaxV_R_minus_L = ONE / ( MaxV_R - MaxV_L );

      for (int v=0; v<NWAVE; v++)
         Flux_Out[ idx_wave[v] ] = _MaxV_R_minus_L*( MaxV_R*Flux_L[ idx_wave[v] ] - MaxV_L*Flux_R[ idx_wave[v] ] );
   }

// longitudinal magnetic flux is always zero
#  ifdef MHD
   Flux_Out[MAG_OFFSET] = ZERO;
#  endif


// 5. evaluate the fluxes of passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out[FLUX_DENS] >= ZERO )
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoL;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
   }

   else
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoR;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
   }
#  endif


// 7. restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_RiemannSolver_HLLE



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLE__
