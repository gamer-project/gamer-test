#ifndef __CUFLU_RIEMANNSOLVER_HLLC__
#define __CUFLU_RIEMANNSOLVER_HLLC__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real ArrayIn[] );

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer extended to include the contact wave
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro (1999)
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput., 18, 1553
//                         c. Coleman, M. S. B. 2020, ApJS, 248, 7
//                3. Wave-speed estimator is set by HLLC_WAVESPEED in CUFLU.h
//                4. Support general EoS
//                5. Shared by the MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out          : Array to store the output flux
//                L/R_In            : Input left/right states (conserved variables)
//                MinDens/Pres      : Density and pressure floors
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_GuessHTilde   :           . . .        the guessed reduced specific enthalpy (SRHD only)
//                EoS_HTilde2Temp   :           . . .        the temperature (SRHD only)
//                EoS_Temper2CSqr   :           . . .        the sound speed square (SRHD only)
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde, 
                               const EoS_H2TEM_t EoS_HTilde2Temp, const EoS_TEM2C_t EoS_Temper2CSqr,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] )
{

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL], R[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );

#  ifdef SRHD
   real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL]; /* primitive vars. */
   real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
   real Usl[NCOMP_TOTAL], Usr[NCOMP_TOTAL];
   real cslsq, csrsq, gammasql, gammasqr;
   real ssl, ssr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
   real lmdal,lmdar; /* Left and Right wave speeds */
   real lmdas; /* Contact wave speed */
   real a,b,c;
   real den,ps; /* Pressure in inner region */
   real lV1, rV1, lV2, rV2, lV3, rV3;
   real lFactor,rFactor; /* Lorentz factor */


/*  1. compute primitive vars. from conserved vars. */
    Hydro_Con2Pri (L, PL, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, 
                  (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, &lFactor );

    Hydro_Con2Pri (R, PR, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                   (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, &rFactor );

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    Hydro_CheckUnphysical(NULL, PL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
    Hydro_CheckUnphysical(NULL, PR, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

/*  2. Transform 4-velocity to 3-velocity */
    lV1=PL[1]/lFactor;
    lV2=PL[2]/lFactor;
    lV3=PL[3]/lFactor;

    rV1=PR[1]/rFactor;
    rV2=PR[2]/rFactor;
    rV3=PR[3]/rFactor;



/*  3. Compute the max and min wave speeds used in Mignone */
    cslsq = EoS_Temper2CSqr( PL[0], PL[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
    csrsq = EoS_Temper2CSqr( PR[0], PR[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    if ( cslsq >= 1.0 || csrsq >= 1.0 || cslsq < 0.0 || csrsq < 0.0 )
      printf( "cslsq=%10.7e, cslrq=%10.7e\n", cslsq, csrsq);
#   endif


//  square of Lorentz factor
    gammasql = SQR(lFactor);
    gammasqr = SQR(rFactor);

    ssl = cslsq / FMA( - gammasql, cslsq, gammasql ); /* Mignone Eq 22.5 */
    ssr = csrsq / FMA( - gammasqr, csrsq, gammasqr ); /* Mignone Eq 22.5 */

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    if ( ( ssl < (real)0.0 ) || ( ssr < (real)0.0 ) ) printf("ssl = %14.7e, ssr = %14.7e\n", ssl, ssr);
#   endif

    real lV2s = lV2*lV2;
    real rV2s = rV2*rV2;

    real lV3s = lV3*lV3;
    real rV3s = rV3*rV3;
 
    real __gammasql = (real)1.0 / gammasql;
    real __gammasqr = (real)1.0 / gammasqr;

    real deltal = ssl*ssl + ssl*( __gammasql + lV2s + lV3s );
    real deltar = ssr*ssr + ssr*( __gammasqr + rV2s + rV3s );

    real ssl__ = (real)1.0 + ssl;
    real ssr__ = (real)1.0 + ssr;


    lmdapl = ( lV1 + SQRT(deltal) ) / ssl__ ;
    lmdaml = ( lV1 - SQRT(deltal) ) / ssl__ ;

    lmdapr = ( rV1 + SQRT(deltar) ) / ssr__ ;
    lmdamr = ( rV1 - SQRT(deltar) ) / ssr__ ;

    lmdal = FMIN(lmdaml, lmdamr); /* Mignone Eq 21 */
    lmdar = FMAX(lmdapl, lmdapr);

/*  4. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18) 
 *     compute HLL conserved quantities using Mignone eq 9
 *  */
    Fl[0] = L[0] * lV1;
    Fl[1] = FMA( L[1], lV1, PL[4] );
    Fl[2] = L[2] * lV1;
    Fl[3] = L[3] * lV1;
#   ifdef REDUCED_ENERGY
    Fl[4] = ( L[4] + PL[4] ) * lV1;
#   else
    Fl[4] = L[1];
#   endif

    if( lmdal >= (real)0.0)
    { /* Fl */
      /* intercell flux is left flux */
      Flux_Out[0] = Fl[0];
      Flux_Out[1] = Fl[1];
      Flux_Out[2] = Fl[2];
      Flux_Out[3] = Fl[3];
      Flux_Out[4] = Fl[4];

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET  );
      return;
    }

    Fr[0] = R[0] * rV1;
    Fr[1] = FMA( R[1], rV1, PR[4] );
    Fr[2] = R[2] * rV1;
    Fr[3] = R[3] * rV1;
#   ifdef REDUCED_ENERGY
    Fr[4] = ( R[4] + PR[4] ) * rV1;
#   else
    Fr[4] = R[1];
#   endif

    if( lmdar <= (real)0.0 )
    { /* Fr */
      /* intercell flux is right flux */
      Flux_Out[0] = Fr[0];
      Flux_Out[1] = Fr[1];
      Flux_Out[2] = Fr[2];
      Flux_Out[3] = Fr[3];
      Flux_Out[4] = Fr[4];

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET  );
      return;
    }


/* 6. Compute contact wave speed using larger root from Mignone Eq 18
 *    Physical root is the root with the minus sign
 */
    lmdatlmda = lmdal*lmdar; 

  /* quadratic formuLa calcuLation */
#   ifdef REDUCED_ENERGY
    a = lmdar * L[1] 
      - lmdal * R[1] 
  	  + lmdatlmda*( R[4] + R[0] - L[4] - L[0] );
#   else
    a = lmdar * ( L[1] )
      - lmdal * ( R[1] )
	  + lmdatlmda*( R[4] - L[4] );
#   endif

#   ifdef REDUCED_ENERGY
    b = lmdal * ( L[4] + L[0] ) - L[1]
      - lmdar * ( R[4] + R[0] ) + R[1]
      + lmdal * ( R[1]*rV1 + PR[4] )
      - lmdar * ( L[1]*lV1 + PL[4] )
      - lmdatlmda*( R[1] - L[1] );
#   else
    b = lmdal * ( L[4] ) - ( L[1] )
      - lmdar * ( R[4] ) + ( R[1] )
      + lmdal * ( R[1]*rV1 + PR[4] )
      - lmdar * ( L[1]*lV1 + PL[4] )
      - lmdatlmda*( R[1] - L[1] );
#   endif

    c = lmdar*R[1] - lmdal*L[1] - ( R[1]*rV1 + PR[4] ) + ( L[1]*lV1 + PL[4] );

    real delta = FMA( b, b, -(real)4*a*c );

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    if (delta < (real) 0.0) printf("delta=%f\n", delta);
#   endif

    lmdas = - ((real)2.0 * c) / ( b + SIGN(b) * SQRT( delta ) );

#   ifdef REDUCED_ENERGY
	ps = lmdas*( ( R[4] + R[0] )*( rV1 - lmdar ) + PR[4]*rV1 ) - R[1]*(rV1 - lmdar) - PR[4];
    ps /= ( lmdas*lmdar - (real)1.0 );
#   else
	ps = lmdas*(  R[4]           *( rV1 - lmdar ) + PR[4]*rV1 ) - R[1]*(rV1 - lmdar) - PR[4];
    ps /= ( lmdas*lmdar - (real)1.0 );
#   endif

	//ps = lmdas*( ( L[4] + L[0] )*( lV1 - lmdal ) + PL[4]*lV1 ) - L[1]*(lV1 - lmdal) - PL[4];
    //ps /= ( lmdas*lmdal - (real)1.0 );


 /* 7. Determine intercell flux according to Mignone 13
 */
   if( lmdas >= (real)0.0 )
   { /* Fls */


    /* now calculate Usl with Mignone Eq 16 */
    den = (real)1.0 / (lmdal - lmdas);

    real factor0 = lmdal - lV1;
    real factor1 = FMA( lmdal, den, -lV1*den );

    Usl[0] =  L[0] * factor1;
    Usl[1] = FMA( L[1], factor0, ps - PL[4] )* den;
    Usl[2] =  L[2] * factor1;
    Usl[3] =  L[3] * factor1;
    Usl[4] = FMA( - PL[4], lV1, FMA( L[4], factor0, ps * lmdas ) ) * den;

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    Hydro_CheckUnphysical(Usl, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdal, Usl[0] - L[0], Fl[0] );
    Flux_Out[1] = FMA( lmdal, Usl[1] - L[1], Fl[1] );
    Flux_Out[2] = FMA( lmdal, Usl[2] - L[2], Fl[2] );
    Flux_Out[3] = FMA( lmdal, Usl[3] - L[3], Fl[3] );
    Flux_Out[4] = FMA( lmdal, Usl[4] - L[4], Fl[4] );

    Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET  );
    return;
   }
   else
   { /* Frs */
    /* now calculate Usr with Mignone Eq 16 */
    den = (real)1.0 / (lmdar - lmdas);

    real factor0 = lmdar - rV1;
    real factor1 = FMA( lmdar, den, -rV1*den );

    Usr[0] = R[0] * factor1;
    Usr[1] = FMA( R[1], factor0, ps - PR[4] ) * den;
    Usr[2] = R[2] * factor1;
    Usr[3] = R[3] * factor1;
    Usr[4] = FMA( - PR[4], rV1, FMA( R[4], factor0, ps * lmdas ) ) * den;

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    Hydro_CheckUnphysical(Usr, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdar, Usr[0] - R[0], + Fr[0] );
    Flux_Out[1] = FMA( lmdar, Usr[1] - R[1], + Fr[1] );
    Flux_Out[2] = FMA( lmdar, Usr[2] - R[2], + Fr[2] );
    Flux_Out[3] = FMA( lmdar, Usr[3] - R[3], + Fr[3] );
    Flux_Out[4] = FMA( lmdar, Usr[4] - R[4], + Fr[4] );

    Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET  );
    return;
   }
#  else //#  ifdef SRHD

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(L[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[0], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(R[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif


// 2. estimate the maximum wave speeds
// 2-1. compute the left/right states
   const real ZERO             = (real)0.0;
   const real ONE              = (real)1.0;
   const real _TWO             = (real)0.5;
   const bool CheckMinPres_Yes = true;
   const real Emag             = NULL_REAL;

   real _RhoL, _RhoR, u_L, u_R, P_L, P_R, Cs_L, Cs_R, W_L, W_R;

   _RhoL = ONE / L[0];
   _RhoR = ONE / R[0];
   u_L   = _RhoL*L[1];
   u_R   = _RhoR*R[1];
   P_L   = Hydro_Con2Pres( L[0], L[1], L[2], L[3], L[4], L+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag,
                           EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   P_R   = Hydro_Con2Pres( R[0], R[1], R[2], R[3], R[4], R+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag,
                           EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   Cs_L  = SQRT(  EoS_DensPres2CSqr( L[0], P_L, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )  );
   Cs_R  = SQRT(  EoS_DensPres2CSqr( R[0], P_R, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )  );

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(P_L) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_L, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(P_R) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_R, __FILE__, __LINE__, __FUNCTION__ );
#  endif


// 2-2a. use the Roe average eigenvalues
#  if   ( HLLC_WAVESPEED == HLL_WAVESPEED_ROE )
   real H_L, H_R, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real u_Roe, v_Roe, w_Roe, H_Roe;

// Roe averages
   H_L             = ( L[4] + P_L )*_RhoL;
   H_R             = ( R[4] + P_R )*_RhoR;
   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   _RhoL_sqrt      = ONE / RhoL_sqrt;
   _RhoR_sqrt      = ONE / RhoR_sqrt;
   _RhoLR_sqrt_sum = ONE / ( RhoL_sqrt + RhoR_sqrt );

   u_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   H_Roe = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

// sound speed
//###NOTE: we have assumed a constant-gamma EoS here
// --> otherwise, one needs to specify how to convert (H-0.5*V2, Rho) to Cs^2
// --> see Eq. [A4] in Coleman 2020
#  if ( EOS != EOS_GAMMA )
#     error : ERROR : HLL_WAVESPEED_ROE only works with EOS_GAMMA !!
#  endif

   const real Gamma    = (real)EoS_AuxArray_Flt[0];
   const real Gamma_m1 = (real)EoS_AuxArray_Flt[1];
   const real _Gamma   = (real)EoS_AuxArray_Flt[3];

   real V2_Roe, Cs2_Roe, Cs_Roe, TempRho, TempPres;

   V2_Roe   = SQR( u_Roe ) + SQR( v_Roe ) + SQR( w_Roe );
   Cs2_Roe  = Gamma_m1*( H_Roe - _TWO*V2_Roe );
   TempRho  = _TWO*( L[0] + R[0] );  // different from Roe average density and is only for applying pressure floor
   TempPres = Cs2_Roe*TempRho*_Gamma;
   TempPres = Hydro_CheckMinPres( TempPres, MinPres );
   Cs2_Roe  = Gamma*TempPres/TempRho;
   Cs_Roe   = SQRT( Cs2_Roe );

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(Cs2_Roe) )
      printf( "ERROR : invalid Cs2_Roe (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Cs2_Roe, __FILE__, __LINE__, __FUNCTION__ );
#  endif

// maximum and minimum eigenvalues
   const real EVal_min = u_Roe - Cs_Roe;
   const real EVal_max = u_Roe + Cs_Roe;

// left/right maximum wave speeds
   W_L = FMIN( EVal_min, u_L-Cs_L );
   W_R = FMAX( EVal_max, u_R+Cs_R );


// 2-2b. use the primitive variable Riemann solver (PVRS)
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_PVRS )
   real Rho_PVRS, Cs_PVRS, RhoCs_PVRS, P_PVRS, Gamma_SL, Gamma_SR, q_L, q_R;

   Rho_PVRS    = _TWO*( L[0] + R[0] );
   Cs_PVRS     = _TWO*( Cs_L + Cs_R );
   RhoCs_PVRS  = Rho_PVRS * Cs_PVRS;
   P_PVRS      = _TWO*(  ( P_L + P_R ) + ( u_L - u_R )*RhoCs_PVRS  );
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
   real u_PVRS, Rho_Cs_PVRS, Rho_SL, Rho_SR, _P;

   u_PVRS      = _TWO*(  ( u_L + u_R ) + ( P_L - P_R )/RhoCs_PVRS  );
   Rho_Cs_PVRS = Rho_PVRS / Cs_PVRS;
   Rho_SL      = L[0] + ( u_L - u_PVRS )*Rho_Cs_PVRS;
   Rho_SR      = R[0] + ( u_PVRS - u_R )*Rho_Cs_PVRS;
   Rho_SL      = FMAX( Rho_SL, MinDens );
   Rho_SR      = FMAX( Rho_SR, MinDens );
   _P          = ONE / P_PVRS;
// see Eq. [9.8] in Toro 1999 for passive scalars
   Gamma_SL    = EoS_DensPres2CSqr( Rho_SL, P_PVRS, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )*Rho_SL*_P;
   Gamma_SR    = EoS_DensPres2CSqr( Rho_SR, P_PVRS, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )*Rho_SR*_P;
#  endif // EOS

   q_L = ( P_PVRS <= P_L ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SL + ONE )/Gamma_SL*( P_PVRS/P_L - ONE )  );
   q_R = ( P_PVRS <= P_R ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SR + ONE )/Gamma_SR*( P_PVRS/P_R - ONE )  );
   W_L = u_L - Cs_L*q_L;
   W_R = u_R + Cs_R*q_R;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(q_L) )
      printf( "ERROR : invalid q_L (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              q_L, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(q_R) )
      printf( "ERROR : invalid q_R (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              q_R, __FILE__, __LINE__, __FUNCTION__ );
#  endif


// 2-2c. use the min/max of the left and right eigenvalues
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_DAVIS )
   const real W_L1 = u_L - Cs_L;
   const real W_L2 = u_R - Cs_R;
   const real W_R1 = u_L + Cs_L;
   const real W_R2 = u_R + Cs_R;
   W_L = FMIN( W_L1, W_L2 );
   W_R = FMAX( W_R1, W_R2 );


#  else
#  error : ERROR : unsupported HLLC_WAVESPEED !!
#  endif // HLLC_WAVESPEED


// 3. evaluate the star-region velocity (V_S) and pressure (P_S)
   real temp1_L, temp1_R, temp2, V_S, P_S;

// do not use u_L-W_L and u_R-W_R to prevent from large round-off errors when Cs<<u~W
// ==> temp1_L ~ temp1_R ~ 0.0 ==> temp2 = inf
#  if   ( HLLC_WAVESPEED == HLL_WAVESPEED_ROE )
   temp1_L = L[0]*(  (EVal_min<u_L-Cs_L) ? (u_L-EVal_min) : (+Cs_L)  );
   temp1_R = R[0]*(  (EVal_max>u_R+Cs_R) ? (u_R-EVal_max) : (-Cs_R)  );
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_PVRS )
   temp1_L = +L[0]*( Cs_L*q_L );
   temp1_R = -R[0]*( Cs_R*q_R );
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_DAVIS )
   temp1_L = +L[0]*(  ( W_L1 < W_L2 ) ? Cs_L : (u_L-u_R)+Cs_R  );
   temp1_R = -R[0]*(  ( W_R2 > W_R1 ) ? Cs_R : (u_L-u_R)+Cs_L  );
#  else
#  error : ERROR : unsupported HLLC_WAVESPEED !!
#  endif

   temp2 = ONE / ( temp1_L - temp1_R );
   V_S   = temp2*( P_L - P_R + temp1_L*u_L - temp1_R*u_R );
   P_S   = temp2*(  temp1_L*( P_R + temp1_R*u_R ) - temp1_R*( P_L + temp1_L*u_L )  );
   P_S   = Hydro_CheckMinPres( P_S, MinPres );


// 4. evaluate the weightings of the left/right fluxes and contact wave
   real Flux_LR[NCOMP_TOTAL], temp4, Coeff_LR, Coeff_S;  // use NCOMP_TOTAL for Flux_LR since it will be passed to Hydro_Con2Flux()

   if ( V_S >= ZERO )
   {
      const real MaxV_L = FMIN( W_L, ZERO );

      Hydro_Con2Flux( 0, Flux_LR, L, MinPres, NULL, NULL, NULL, NULL, &P_L );

      for (int v=0; v<NCOMP_FLUID; v++)   Flux_LR[v] -= MaxV_L*L[v];    // fluxes along the maximum wave speed

//    deal with the special case of V_S=MaxV_L=0
//###REVISE: should it return zero flux due to symmetry?
      if ( V_S == ZERO  &&  MaxV_L == ZERO )
      {
         Coeff_LR = ONE;
         Coeff_S  = ZERO;
      }

      else
      {
         temp4    = ONE / ( V_S - MaxV_L );
         Coeff_LR = temp4*V_S;
         Coeff_S  = -temp4*MaxV_L*P_S;
      }
   } // if ( V_S >= ZERO )

   else // V_S < 0.0
   {
      const real MaxV_R = FMAX( W_R, ZERO );

      Hydro_Con2Flux( 0, Flux_LR, R, MinPres, NULL, NULL, NULL, NULL, &P_R );

      for (int v=0; v<NCOMP_FLUID; v++)    Flux_LR[v] -= MaxV_R*R[v];   // fluxes along the maximum wave speed

      temp4    = ONE / ( V_S - MaxV_R );
      Coeff_LR = temp4*V_S;
      Coeff_S  = -temp4*MaxV_R*P_S;
   } // if ( V_S >= ZERO ) ... else ...


// 5. evaluate the HLLC fluxes
   for (int v=0; v<NCOMP_FLUID; v++)   Flux_Out[v] = Coeff_LR*Flux_LR[v];

   Flux_Out[1] += Coeff_S;
   Flux_Out[4] += Coeff_S*V_S;


// 6. evaluate the fluxes of passive scalars
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

#  endif // #  ifdef SRHD



} // FUNCTION : Hydro_RiemannSolver_HLLC



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLC__
