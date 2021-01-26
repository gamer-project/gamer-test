#ifndef __CUFLU_RIEMANNSOLVER_ROE__
#define __CUFLU_RIEMANNSOLVER_ROE__



#include "CUFLU.h"

#if ( MODEL == HYDRO &&  !defined SRHD)



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#if   ( CHECK_INTERMEDIATE == EXACT )
# include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( CHECK_INTERMEDIATE == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( CHECK_INTERMEDIATE == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#elif ( CHECK_INTERMEDIATE == HLLD )
# include "CUFLU_Shared_RiemannSolver_HLLD.cu"
#endif

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real ArrayIn[] );
#if   ( CHECK_INTERMEDIATE == EXACT )
void Hydro_RiemannSolver_Exact( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                const real Dens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                                const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#elif ( CHECK_INTERMEDIATE == HLLE )
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Dens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#elif ( CHECK_INTERMEDIATE == HLLC )
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp, const EoS_TEM2C_t EoS_Temper2CSqr,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                               const real* const EoS_Table[EOS_NTABLE_MAX] );
#elif ( CHECK_INTERMEDIATE == HLLD )
void Hydro_RiemannSolver_HLLD( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Dens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_Roe
// Description :  Approximate Riemann solver of Roe
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : (a) "Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro"
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                3. Shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out          : Array to store the output flux
//                L_In              : Input left  state (conserved variables)
//                R_In              : Input right state (conserved variables)
//                MinDens/Pres      : Density and pressure floors
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                              const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                              const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                              const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] )
{

// check
#  if ( EOS == EOS_GAMMA )
   const real *Passive = NULL;   // EOS_GAMMA does not involve passive scalars
#  elif ( defined GAMER_DEBUG )
   printf( "ERROR : EOS != EOS_GAMMA is NOT supported at file <%s>, line <%d>, function <%s> !!\n",
           __FILE__, __LINE__, __FUNCTION__ );
#  endif


// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL_PLUS_MAG], R[NCOMP_TOTAL_PLUS_MAG];

   for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );

// longitudinal B field in the left and right states should be the same
#  if ( defined GAMER_DEBUG  &&  defined MHD )
   if ( L[MAG_OFFSET] != R[MAG_OFFSET] )
      printf( "ERROR : BxL (%24.17e) != BxR (%24.17e) for XYZ %d at file <%s>, line <%d>, function <%s>!!\n",
              L[MAG_OFFSET], R[MAG_OFFSET], XYZ, __FILE__, __LINE__, __FUNCTION__ );
#  endif


// 2. evaluate the average values
   const real ZERO             = (real)0.0;
   const real ONE              = (real)1.0;
   const real _TWO             = (real)0.5;
   const real Gamma            = EoS_AuxArray_Flt[0];    // only support constant-gamma EoS (i.e., EOS_GAMMA)
   const real Gamma_m1         = EoS_AuxArray_Flt[1];
   const bool CheckMinPres_Yes = true;
#  ifdef MHD
   const real TWO              = (real)2.0;
   const real Gamma_m2         = Gamma - TWO;
#  endif

   real Rho, _Rho, _RhoL, _RhoR, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real PL, PR, HL, HR, u, v, w, V2, H, a, a2, GammaP_Rho, EmagL=NULL_REAL, EmagR=NULL_REAL;
#  ifdef MHD
   real Rho_sqrt, _Rho_sqrt;                    // Roe-average density
   real ByL, BzL, ByR, BzR;                     // magnetic field from left and right states
   real Bx, By, Bz, B2, Bn2, Bn, B2_Rho;        // Roe-average magnetic field
   real Cax, Cax2, Cat2, Cs, Cs2, Cf, Cf2;      // Alfven, slow, and fast waves
   real alpha_f, alpha_s, beta_y, beta_z;       // Eqs. (A16) and (A17) in ref-b
   real Ca2_plus_a2, Ca2_min_a2, Cf2_min_Cs2;   // Ca^2+a^2, Ca^2-a^2, Cf^2-Cs^2
   real X, Y;                                   // Eqs. (B15) and (B16) in ref-b
   real S;                                      // sign(Bx)
   real Bn_star;                                // Eq. (B20) in ref-b
   real beta_n_star2, beta_y_star, beta_z_star; // Eq. (B28) in ref-b
#  endif

   _RhoL = ONE / L[0];
   _RhoR = ONE / R[0];
#  ifdef MHD
   Bx    = L[ MAG_OFFSET + 0 ];  // assuming Bx=BxL=BxR
   ByL   = L[ MAG_OFFSET + 1 ];
   BzL   = L[ MAG_OFFSET + 2 ];
   ByR   = R[ MAG_OFFSET + 1 ];
   BzR   = R[ MAG_OFFSET + 2 ];
   EmagL = _TWO*( SQR(Bx) + SQR(ByL) + SQR(BzL) );
   EmagR = _TWO*( SQR(Bx) + SQR(ByR) + SQR(BzR) );
#  endif
   PL    = Hydro_Con2Pres( L[0], L[1], L[2], L[3], L[4], L+NCOMP_FLUID, CheckMinPres_Yes, MinPres, EmagL,
                           EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   PR    = Hydro_Con2Pres( R[0], R[1], R[2], R[3], R[4], R+NCOMP_FLUID, CheckMinPres_Yes, MinPres, EmagR,
                           EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
#  ifdef MHD
   HL    = _RhoL*( L[4] + PL + EmagL );
   HR    = _RhoR*( R[4] + PR + EmagR );
#  else
   HL    = _RhoL*( L[4] + PL );
   HR    = _RhoR*( R[4] + PR );
#  endif

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(L[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[0], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(R[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   Rho             = RhoL_sqrt*RhoR_sqrt;
   _Rho            = ONE/Rho;
   _RhoL_sqrt      = ONE/RhoL_sqrt;
   _RhoR_sqrt      = ONE/RhoR_sqrt;
   _RhoLR_sqrt_sum = ONE/(RhoL_sqrt + RhoR_sqrt);

   u  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   V2 = u*u + v*v + w*w;
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*HL   +  RhoR_sqrt*HR   );

#  ifdef MHD
   Rho_sqrt  = SQRT( Rho );
   _Rho_sqrt = ONE/Rho_sqrt;
   By        = _RhoLR_sqrt_sum*( RhoL_sqrt*ByR + RhoR_sqrt*ByL );
   Bz        = _RhoLR_sqrt_sum*( RhoL_sqrt*BzR + RhoR_sqrt*BzL );
   Bn2       = SQR( By ) + SQR( Bz );
   Bn        = SQRT( Bn2 );
   B2        = SQR( Bx ) + Bn2;
   B2_Rho    = B2*_Rho;
   S         = SIGN( Bx );
   X         = _TWO*( SQR(ByR-ByL) + SQR(BzR-BzL) )*SQR( _RhoLR_sqrt_sum );
   X        *= Gamma_m2;
#  ifdef EULERY
   Y         = _TWO*( L[0] + R[0] )*_Rho ;
#  else
   Y         = ONE;
#  endif
   Y        *= Gamma_m2;
   Bn_star   = SQRT( Gamma_m1 - Y )*Bn;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(Gamma_m1-Y) )
      printf( "ERROR : invalid Gamma_m1-Y (%14.7e, Gamma_m1 %14.7e, Y %14.7e) at file <%s>, line <%d>, function <%s>\n",
              Gamma_m1-Y, Gamma_m1, Y, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   if ( Bn == ZERO ) {
      beta_y = ONE;
      beta_z = ZERO;
   }
   else {
      const real _Bn = ONE/Bn;
      beta_y = By*_Bn;
      beta_z = Bz*_Bn;
   }

   if ( Bn_star == ZERO ) {
      beta_y_star = ONE;
      beta_z_star = ZERO;
   }
   else {
      const real _Bn_star = ONE/Bn_star;
      beta_y_star = By*_Bn_star;
      beta_z_star = Bz*_Bn_star;
   }

   beta_n_star2 = SQR( beta_y_star ) + SQR( beta_z_star );
#  endif // #ifdef MHD

   GammaP_Rho = Gamma_m1*( H - _TWO*V2 );
#  ifdef MHD
   GammaP_Rho -= Gamma_m1*B2_Rho;   // H = 0.5*v^2 + B^2/rho + gamma/(gamma-1)*P/rho
#  endif
   GammaP_Rho = Gamma*_Rho*Hydro_CheckMinPres( GammaP_Rho*Rho/Gamma, MinPres );  // apply pressure floor

   a2  = GammaP_Rho;
#  ifdef MHD
   a2 -= X;
#  endif
#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(a2) )
      printf( "ERROR : invalid a2 (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              a2, __FILE__, __LINE__, __FUNCTION__ );
#  endif
   a  = SQRT( a2 );

#  ifdef MHD
   Cax2        = SQR(Bx)*_Rho;
   Cax         = SQRT( Cax2 );
   Cat2        = ( Gamma_m1 - Y )*Bn2*_Rho;
   Ca2_plus_a2 = Cat2 + Cax2 + a2;
   Ca2_min_a2  = Cat2 + Cax2 - a2;
   Cf2_min_Cs2 = SQRT( SQR(Ca2_min_a2) + (real)4.0*a2*Cat2 );

// evaluate the fast/slow wave speed (Cf/Cs)
   if ( Cat2 == ZERO )
   {
      if ( Cax2 == a2 ) {
         Cf2 = a2;
         Cs2 = a2;
      }
      else if ( Cax2 > a2 ) {
         Cf2 = Cax2;
         Cs2 = a2;
      }
      else {
         Cf2 = a2;
         Cs2 = Cax2;
      }
   }

   else
   {
      if ( Cax2 == ZERO ) {
         Cf2 = a2 + Cat2;
         Cs2 = ZERO;
      }
      else {
         Cf2 = _TWO*( Ca2_plus_a2 + Cf2_min_Cs2 );
         Cs2 = a2*Cax2/Cf2;   // do not use "Cf2 - Cf2_min_Cs2" to avoid negative values caused by round-off errors
//       Cs2 = Cf2 - Cf2_min_Cs2;
      }
   } // if ( Cat2 == ZERO ) ... else ...

   Cf = SQRT( Cf2 );
   Cs = SQRT( Cs2 );

   const real a2_min_Cs2 = a2 - Cs2;
   const real Cf2_min_a2 = Cf2 - a2;

   if ( Cf2_min_Cs2 == ZERO ) {
      alpha_f = ONE;
      alpha_s = ZERO;
   }
   else if ( a2_min_Cs2 <= ZERO ) {
      alpha_f = ZERO;
      alpha_s = ONE;
   }
   else if ( Cf2_min_a2 <= ZERO ) {
      alpha_f = ONE;
      alpha_s = ZERO;
   }
   else {
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(a2_min_Cs2) )
         printf( "ERROR : invalid a2_min_Cs2 (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 a2_min_Cs2, __FILE__, __LINE__, __FUNCTION__ );

      if ( Hydro_CheckNegative(Cf2_min_a2) )
         printf( "ERROR : invalid Cf2_min_a2 (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Cf2_min_a2, __FILE__, __LINE__, __FUNCTION__ );
#     endif

      const real _Cf2_min_Cs2 = ONE/Cf2_min_Cs2;
      alpha_f = SQRT( a2_min_Cs2*_Cf2_min_Cs2 );
      alpha_s = SQRT( Cf2_min_a2*_Cf2_min_Cs2 );
   }
#  endif // #ifdef MHD


// 3. evaluate the eigenvalues
#  ifdef MHD
   const real EigenVal[NWAVE] = { u-Cf, u-Cax, u-Cs, u, u+Cs, u+Cax, u+Cf };
#  else
   const real EigenVal[NWAVE] = { u-a, u, u, u, u+a };
#  endif


// 4. evaluate the left and right fluxes
   real Flux_L[NCOMP_TOTAL_PLUS_MAG], Flux_R[NCOMP_TOTAL_PLUS_MAG];

   Hydro_Con2Flux( 0, Flux_L, L, MinPres, NULL, NULL, NULL, NULL, &PL );
   Hydro_Con2Flux( 0, Flux_R, R, MinPres, NULL, NULL, NULL, NULL, &PR );

// 5. return the upwind fluxes if flow is supersonic
   if ( EigenVal[0] >= ZERO )
   {
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   Flux_Out[v] = Flux_L[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

      return;
   }

   if ( EigenVal[NWAVE-1] <= ZERO )
   {
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   Flux_Out[v] = Flux_R[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

      return;
   }


// 6. evaluate the eigenvectors
//    --> right eigenvectors: columns of REigenVec[][]
//        left  eigenvectors: rows    of LEigenVec[][]
#  ifdef MHD
   real REigenVec[NWAVE][NWAVE], LEigenVec[NWAVE][NWAVE];

// Eqs. [A13]-[A17] in ref-b
   const real Af  = a*alpha_f*Rho_sqrt;
   const real As  = a*alpha_s*Rho_sqrt;
   const real Cff = Cf*alpha_f;
   const real Css = Cs*alpha_s;
   const real Qf  = Cff*S;
   const real Qs  = Css*S;

// right eigenvectors
   REigenVec[0][0] = alpha_f;
   REigenVec[0][1] = ZERO;
   REigenVec[0][2] = alpha_s;
   REigenVec[0][3] = ONE;
   REigenVec[0][4] = alpha_s;
   REigenVec[0][5] = ZERO;
   REigenVec[0][6] = alpha_f;

   const real u_alpha_f = u*alpha_f;
   const real u_alpha_s = u*alpha_s;

   REigenVec[1][0] = u_alpha_f - Cff;
   REigenVec[1][1] = ZERO;
   REigenVec[1][2] = u_alpha_s - Css;
   REigenVec[1][3] = u;
   REigenVec[1][4] = u_alpha_s + Css;
   REigenVec[1][5] = ZERO;
   REigenVec[1][6] = u_alpha_f + Cff;

   const real v_alpha_f      = v*alpha_f;
   const real v_alpha_s      = v*alpha_s;
   const real Qs_beta_y_star = Qs*beta_y_star;
   const real Qf_beta_y_star = Qf*beta_y_star;

   REigenVec[2][0] = v_alpha_f + Qs_beta_y_star;
   REigenVec[2][1] = -beta_z;
   REigenVec[2][2] = v_alpha_s - Qf_beta_y_star;
   REigenVec[2][3] = v;
   REigenVec[2][4] = v_alpha_s + Qf_beta_y_star;
   REigenVec[2][5] = beta_z;
   REigenVec[2][6] = v_alpha_f - Qs_beta_y_star;

   const real w_alpha_f      = w*alpha_f;
   const real w_alpha_s      = w*alpha_s;
   const real Qs_beta_z_star = Qs*beta_z_star;
   const real Qf_beta_z_star = Qf*beta_z_star;

   REigenVec[3][0] = w_alpha_f + Qs_beta_z_star;
   REigenVec[3][1] = beta_y;
   REigenVec[3][2] = w_alpha_s - Qf_beta_z_star;
   REigenVec[3][3] = w;
   REigenVec[3][4] = w_alpha_s + Qf_beta_z_star;
   REigenVec[3][5] = -beta_y;
   REigenVec[3][6] = w_alpha_f - Qs_beta_z_star;

   const real Hp           = H - B2_Rho;  // H_prime = H - B^2/rho
   const real u_Cf         = u*Cf;
   const real u_Cs         = u*Cs;
   const real v_by_w_bz    = v*beta_y_star + w*beta_z_star;
   const real Qs_v_by_w_bz = Qs*v_by_w_bz;
   const real Qf_v_by_w_bz = Qf*v_by_w_bz;
   const real Bn_b2_Rho    = Bn_star*beta_n_star2*_Rho;
   const real As_Bn_b2_Rho = As*Bn_b2_Rho;
   const real Af_Bn_b2_Rho = Af*Bn_b2_Rho;

   REigenVec[4][0] = alpha_f*( Hp - u_Cf ) + Qs_v_by_w_bz + As_Bn_b2_Rho;
   REigenVec[4][1] = -( v*beta_z - w*beta_y );
   REigenVec[4][2] = alpha_s*( Hp - u_Cs ) - Qf_v_by_w_bz - Af_Bn_b2_Rho;
   REigenVec[4][3] = _TWO*V2 + X/Gamma_m1;
   REigenVec[4][4] = alpha_s*( Hp + u_Cs ) + Qf_v_by_w_bz - Af_Bn_b2_Rho;
   REigenVec[4][5] = -REigenVec[4][1];
   REigenVec[4][6] = alpha_f*( Hp + u_Cf ) - Qs_v_by_w_bz + As_Bn_b2_Rho;

   const real beta_y_star_Rho = beta_y_star*_Rho;

   REigenVec[5][0] = As*beta_y_star_Rho;
   REigenVec[5][1] = -S*beta_z*_Rho_sqrt;
   REigenVec[5][2] = -Af*beta_y_star_Rho;
   REigenVec[5][3] = ZERO;
   REigenVec[5][4] = REigenVec[5][2];
   REigenVec[5][5] = REigenVec[5][1];
   REigenVec[5][6] = REigenVec[5][0];

   const real beta_z_star_Rho = beta_z_star*_Rho;

   REigenVec[6][0] = As*beta_z_star_Rho;
   REigenVec[6][1] = S*beta_y*_Rho_sqrt;
   REigenVec[6][2] = -Af*beta_z_star_Rho;
   REigenVec[6][3] = ZERO;
   REigenVec[6][4] = REigenVec[6][2];
   REigenVec[6][5] = REigenVec[6][1];
   REigenVec[6][6] = REigenVec[6][0];


// left eigenvectors
   const real Qy_star     = beta_y_star / beta_n_star2;  // Eq. [B30] in ref-b
   const real Qz_star     = beta_z_star / beta_n_star2;
   const real norm_hat    = _TWO / a2;
   const real norm_bar    = norm_hat*Gamma_m1;
   const real Cff_hat     = norm_hat*Cff;
   const real Css_hat     = norm_hat*Css;
   const real Af_hat      = norm_hat*Af;
   const real As_hat      = norm_hat*As;
   const real Qf_hat      = norm_hat*Qf;
   const real Qs_hat      = norm_hat*Qs;
   const real X_hat       = norm_hat*X;
   const real alpha_f_bar = norm_bar*alpha_f;
   const real alpha_s_bar = norm_bar*alpha_s;

   const real alpha_f_V2_Hp  = alpha_f_bar*( V2 - Hp);
   const real alpha_s_V2_Hp  = alpha_s_bar*( V2 - Hp);
   const real v_Qy_w_Qz      = v*Qy_star + w*Qz_star;
   const real Qf_v_Qy_w_Qz   = Qf_hat*v_Qy_w_Qz;
   const real Qs_v_Qy_w_Qz   = Qs_hat*v_Qy_w_Qz;
   const real Bn_star_Rho    = Bn_star*_Rho;
   const real Af_Bn_star_Rho = Af_hat*Bn_star_Rho;
   const real As_Bn_star_Rho = As_hat*Bn_star_Rho;

   LEigenVec[0][0] = alpha_f_V2_Hp + Cff_hat*( Cf + u ) - Qs_v_Qy_w_Qz - As_Bn_star_Rho;
   LEigenVec[1][0] = _TWO*( v*beta_z - w*beta_y );
   LEigenVec[2][0] = alpha_s_V2_Hp + Css_hat*( Cs + u ) + Qf_v_Qy_w_Qz + Af_Bn_star_Rho;
   LEigenVec[3][0] = ONE - norm_bar*V2 + TWO*X_hat;
   LEigenVec[4][0] = alpha_s_V2_Hp + Css_hat*( Cs - u ) - Qf_v_Qy_w_Qz + Af_Bn_star_Rho;
   LEigenVec[5][0] = -LEigenVec[1][0];
   LEigenVec[6][0] = alpha_f_V2_Hp + Cff_hat*( Cf - u ) + Qs_v_Qy_w_Qz - As_Bn_star_Rho;

   const real u_alpha_f_bar = u*alpha_f_bar;
   const real u_alpha_s_bar = u*alpha_s_bar;

   LEigenVec[0][1] = -u_alpha_f_bar - Cff_hat;
   LEigenVec[1][1] = ZERO;
   LEigenVec[2][1] = -u_alpha_s_bar - Css_hat;
   LEigenVec[3][1] = TWO*norm_bar*u;
   LEigenVec[4][1] = -u_alpha_s_bar + Css_hat;
   LEigenVec[5][1] = ZERO;
   LEigenVec[6][1] = -u_alpha_f_bar + Cff_hat;

   const real v_alpha_f_bar = v*alpha_f_bar;
   const real v_alpha_s_bar = v*alpha_s_bar;
   const real Qs_Qy_star    = Qs_hat*Qy_star;
   const real Qf_Qy_star    = Qf_hat*Qy_star;

   LEigenVec[0][2] = -v_alpha_f_bar + Qs_Qy_star;
   LEigenVec[1][2] = -_TWO*beta_z;
   LEigenVec[2][2] = -v_alpha_s_bar - Qf_Qy_star;
   LEigenVec[3][2] = TWO*norm_bar*v;
   LEigenVec[4][2] = -v_alpha_s_bar + Qf_Qy_star;
   LEigenVec[5][2] = -LEigenVec[1][2];
   LEigenVec[6][2] = -v_alpha_f_bar - Qs_Qy_star;

   const real w_alpha_f_bar = w*alpha_f_bar;
   const real w_alpha_s_bar = w*alpha_s_bar;
   const real Qs_Qz_star    = Qs_hat*Qz_star;
   const real Qf_Qz_star    = Qf_hat*Qz_star;

   LEigenVec[0][3] = -w_alpha_f_bar + Qs_Qz_star;
   LEigenVec[1][3] =  _TWO*beta_y;
   LEigenVec[2][3] = -w_alpha_s_bar - Qf_Qz_star;
   LEigenVec[3][3] = TWO*norm_bar*w;
   LEigenVec[4][3] = -w_alpha_s_bar + Qf_Qz_star;
   LEigenVec[5][3] = -LEigenVec[1][3];
   LEigenVec[6][3] = -w_alpha_f_bar - Qs_Qz_star;

   LEigenVec[0][4] = alpha_f_bar;
   LEigenVec[1][4] = ZERO;
   LEigenVec[2][4] = alpha_s_bar;
   LEigenVec[3][4] = -Gamma_m1 / a2;
   LEigenVec[4][4] = alpha_s_bar;
   LEigenVec[5][4] = ZERO;
   LEigenVec[6][4] = alpha_f_bar;

   LEigenVec[0][5] = As_hat*Qy_star - alpha_f_bar*By;
   LEigenVec[1][5] = -_TWO*S*beta_z*Rho_sqrt;
   LEigenVec[2][5] = -Af_hat*Qy_star - alpha_s_bar*By;
   LEigenVec[3][5] = TWO*norm_bar*By;
   LEigenVec[4][5] = LEigenVec[2][5];
   LEigenVec[5][5] = LEigenVec[1][5];
   LEigenVec[6][5] = LEigenVec[0][5];

   LEigenVec[0][6] = As_hat*Qz_star - alpha_f_bar*Bz;
   LEigenVec[1][6] = _TWO*S*beta_y*Rho_sqrt;
   LEigenVec[2][6] = -Af_hat*Qz_star - alpha_s_bar*Bz;
   LEigenVec[3][6] = TWO*norm_bar*Bz;
   LEigenVec[4][6] = LEigenVec[2][6];
   LEigenVec[5][6] = LEigenVec[1][6];
   LEigenVec[6][6] = LEigenVec[0][6];


#  else // #ifdef MHD
   const real REigenVec[NWAVE][NWAVE] =
      {  {   ONE,     ONE, ZERO, ZERO,   ONE },
         {   u-a,       u, ZERO, ZERO,   u+a },
         {     v,       v,  ONE, ZERO,     v },
         {     w,       w, ZERO,  ONE,     w },
         { H-u*a, _TWO*V2,    v,    w, H+u*a }  };
#  endif // #ifdef MHD ... else ...


// 7. evaluate the amplitudes along different characteristics (eigenvectors)

// index mapping between arrays with size NWAVE and NCOMP_TOTAL_PLUS_MAG;
#  ifdef MHD
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4, MAG_OFFSET+1, MAG_OFFSET+2 };
#  else
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };
#  endif
   real Jump[NWAVE], Amp[NWAVE];

   for (int v=0; v<NWAVE; v++)  Jump[v] = R[ idx_wave[v] ] - L[ idx_wave[v] ];

#  ifdef MHD
   Amp[0] =  LEigenVec[0][0]*Jump[0] + LEigenVec[0][1]*Jump[1] + LEigenVec[0][2]*Jump[2] + LEigenVec[0][3]*Jump[3];
   Amp[1] =  LEigenVec[1][0]*Jump[0]                           + LEigenVec[1][2]*Jump[2] + LEigenVec[1][3]*Jump[3];
   Amp[2] =  LEigenVec[2][0]*Jump[0] + LEigenVec[2][1]*Jump[1] + LEigenVec[2][2]*Jump[2] + LEigenVec[2][3]*Jump[3];
   Amp[3] =  LEigenVec[3][0]*Jump[0] + LEigenVec[3][1]*Jump[1] + LEigenVec[3][2]*Jump[2] + LEigenVec[3][3]*Jump[3]
           + LEigenVec[3][4]*Jump[4] + LEigenVec[3][5]*Jump[5] + LEigenVec[3][6]*Jump[6];
   Amp[4] =  LEigenVec[4][0]*Jump[0] + LEigenVec[4][1]*Jump[1] + LEigenVec[4][2]*Jump[2] + LEigenVec[4][3]*Jump[3];
   Amp[5] = -Amp[1];
   Amp[6] =  LEigenVec[6][0]*Jump[0] + LEigenVec[6][1]*Jump[1] + LEigenVec[6][2]*Jump[2] + LEigenVec[6][3]*Jump[3];

   const real tmp0 = LEigenVec[0][4]*Jump[4] + LEigenVec[0][5]*Jump[5] + LEigenVec[0][6]*Jump[6];
   const real tmp1 =                           LEigenVec[1][5]*Jump[5] + LEigenVec[1][6]*Jump[6];
   const real tmp2 = LEigenVec[2][4]*Jump[4] + LEigenVec[2][5]*Jump[5] + LEigenVec[2][6]*Jump[6];

   Amp[0] += tmp0;
   Amp[1] += tmp1;
   Amp[2] += tmp2;
   Amp[4] += tmp2;
   Amp[5] += tmp1;
   Amp[6] += tmp0;

#  else // #ifdef MHD
   Amp[2] = Jump[2] - v*Jump[0];
   Amp[3] = Jump[3] - w*Jump[0];
   Amp[1] = Gamma_m1/a2*( Jump[0]*(H-SQR(u)) + u*Jump[1] - Jump[4] + v*Amp[2] + w*Amp[3] );
   Amp[0] = _TWO/a*( Jump[0]*(u+a) - Jump[1] - a*Amp[1] );
   Amp[4] = Jump[0] - Amp[0] - Amp[1];
#  endif // #ifdef MHD ... else ...


// 8. verify that the density and pressure in the intermediate states are positive
#  ifdef CHECK_INTERMEDIATE
   const bool CheckMinPres_No = false;
   real I_Pres, I_States[ NCOMP_FLUID + NCOMP_MAG ];

   for (int v=0; v<NCOMP_FLUID; v++)   I_States[ v               ] = L[ v              ];
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG;   v++)   I_States[ v + NCOMP_FLUID ] = L[ v + MAG_OFFSET ];
#  endif

   for (int t=0; t<NWAVE-1; t++)
   {
      for (int v=0; v<NCOMP_FLUID; v++)       I_States[ v     ] += Amp[t]*REigenVec[v][t];
#     ifdef MHD
      for (int v=NCOMP_FLUID; v<NWAVE; v++)   I_States[ v + 1 ] += Amp[t]*REigenVec[v][t];
#     endif

      if ( EigenVal[t+1] > EigenVal[t] )  // skip the degenerate states
      {
#        ifdef MHD
         const real Emag = _TWO*(  SQR( I_States[NCOMP_FLUID+0] )
                                 + SQR( I_States[NCOMP_FLUID+1] )
                                 + SQR( I_States[NCOMP_FLUID+2] )  );
#        else
         const real Emag = NULL_REAL;
#        endif
         I_Pres = Hydro_Con2Pres( I_States[0], I_States[1], I_States[2], I_States[3], I_States[4], Passive,
                                  CheckMinPres_No, NULL_REAL, Emag, EoS_DensEint2Pres, NULL, NULL,
                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

//       if unphysical results occur, recalculate fluxes by a substitute Riemann solver
         if ( I_States[0] <= ZERO  ||  I_Pres <= ZERO )
         {
#           ifdef GAMER_DEBUG
            printf( "WARNING : intermediate states check failed (density %14.7e, pressure %14.7e) !!\n",
                    I_States[0], I_Pres );
#           endif

#           if   ( CHECK_INTERMEDIATE == EXACT  &&  !defined MHD )
            Hydro_RiemannSolver_Exact( 0, Flux_Out, L, R, MinDens, MinPres, EoS_DensEint2Pres, EoS_DensPres2CSqr,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#           elif ( CHECK_INTERMEDIATE == HLLE )
            Hydro_RiemannSolver_HLLE ( 0, Flux_Out, L, R, MinDens, MinPres, EoS_DensEint2Pres, EoS_DensPres2CSqr,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#           elif ( CHECK_INTERMEDIATE == HLLC  &&  !defined MHD )
            Hydro_RiemannSolver_HLLC ( 0, Flux_Out, L, R, MinDens, MinPres, EoS_DensEint2Pres, EoS_DensPres2CSqr,
                                       NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#           elif ( CHECK_INTERMEDIATE == HLLD  &&  defined MHD )
            Hydro_RiemannSolver_HLLD ( 0, Flux_Out, L, R, MinDens, MinPres, EoS_DensEint2Pres, EoS_DensPres2CSqr,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

#           else
#           error : ERROR : unsupported CHECK_INTERMEDIATE (EXACT/HLLE/HLLC/HLLD) !!
#           endif // CHECK_INTERMEDIATE

            Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );
            return;

         } // if ( I_States[0] <= ZERO  ||  I_Pres <= ZERO )
      } // if ( EigenVal[t+1] > EigenVal[t] )
   } // for (int t=0; t<NWAVE-1; t++)
#  endif // #ifdef CHECK_INTERMEDIATE


// 9. evaluate the Roe fluxes
   for (int v=0; v<NWAVE; v++)   Amp[v] *= FABS( EigenVal[v] );

   for (int v=0; v<NWAVE; v++)
   {
      const int vv = idx_wave[v];

      Flux_Out[vv] = Flux_L[vv] + Flux_R[vv];

      for (int t=0; t<NWAVE; t++)   Flux_Out[vv] -= Amp[t]*REigenVec[v][t];

      Flux_Out[vv] *= _TWO;
   }

// longitudinal magnetic flux is always zero
#  ifdef MHD
   Flux_Out[MAG_OFFSET] = ZERO;
#  endif


// 10. evaluate the fluxes for passive scalars
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


// 11. restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_RiemannSolver_Roe



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_ROE__
