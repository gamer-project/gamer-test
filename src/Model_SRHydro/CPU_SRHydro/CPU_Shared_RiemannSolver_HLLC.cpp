#include "GAMER.h"
#include "CUFLU.h"
#include "../../../include/CPU_prototypes.h"

#if ( MODEL == SR_HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                The wave speed is estimated by the same formula in HLLE solver
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput.,
//                            18, 1553
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out : Array to store the output flux
//                L_In     : Input left  state (conserved variables)
//                R_In     : Input right state (conserved variables)
//                Gamma    : Ratio of specific heats
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannSolver_HLLC( const int XYZ,
                             real Flux_Out[],
                             const real L_In[],
                             const real R_In[],
                             const real Gamma,
                             const real MinPres )
{
  real CL[NCOMP_TOTAL], CR[NCOMP_TOTAL]; /* conserved vars. */
  real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL]; /* primitive vars. */
  real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
  real Fhll[NCOMP_TOTAL], Uhll[NCOMP_TOTAL];
  real Usl[NCOMP_TOTAL], Usr[NCOMP_TOTAL];
  
  const real Gamma_m1 = Gamma - (real)1.0;
  double rhl, rhr, csl, csr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  double ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  double lmdal,lmdar; /* Left and Right wave speeds */
  double lmdas; /* Contact wave speed */
  double ovlrmll;
  double a,b,c,quad;
  double den,ps; /* Pressure in inner region */
  double lV1, lV2, lV3, rV1, rV2, rV3;
  double lFactor,rFactor; /* Lorentz factor */


/* 0. reorder the input conserved variables for different spatial directions */
   for(int v=0;v<NCOMP_TOTAL;v++){
       CL[v]=L_In[v];
       CR[v]=R_In[v];
   }

   CPU_Rotate3D( CL, XYZ, true );
   CPU_Rotate3D( CR, XYZ, true );

/* 1. compute primitive vars. from conserved vars. */
   CPU_Con2Pri (CL, PL, Gamma);
   CPU_Con2Pri (CR, PR, Gamma);

/* 2. Transform 4-velocity to 3-velocity */
   lFactor=1/SQRT(1+SQR(PL[1])+SQR(PL[2])+SQR(PL[3]));
   rFactor=1/SQRT(1+SQR(PR[1])+SQR(PR[2])+SQR(PR[3]));

   lV1=PL[1]*lFactor;
   lV2=PL[2]*lFactor;
   lV3=PL[3]*lFactor;

   rV1=PR[1]*rFactor;
   rV2=PR[2]*rFactor;
   rV3=PR[3]*rFactor;

/* 3. Compute the max and min wave speeds used in Mignone */
   rhl = PL[0] + PL[4] * Gamma / Gamma_m1; /* Mignone Eq 3.5 */
   rhr = PR[0] + PR[4] * Gamma / Gamma_m1;

   csl = SQRT(Gamma * PL[4] / rhl); /* Mignone Eq 4 */
   csr = SQRT(Gamma * PR[4] / rhr);

   cslsq = SQR(csl);
   csrsq = SQR(csr);

   vsql = SQR(lV1) + SQR(lV2) + SQR(lV3);
   vsqr = SQR(rV1) + SQR(rV2) + SQR(rV3);

   gammasql = 1.0 / (1.0 - vsql);
   gammasqr = 1.0 / (1.0 - vsqr);

   ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
   ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

   radl = SQRT( ssl*(1.0-SQR(lV1)+ssl) ); /* Mignone Eq 23 (radical part) */
   radr = SQRT( ssr*(1.0-SQR(rV1)+ssr) );

   lmdapl = (lV1 + radl) / (1.0 + ssl); /* Mignone Eq 23 */
   lmdapr = (rV1 + radr) / (1.0 + ssr);
   lmdaml = (lV1 - radl) / (1.0 + ssl);
   lmdamr = (rV1 - radr) / (1.0 + ssr);

   lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
   lmdar = MAX(lmdapl, lmdapr);
    
/* 4. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18) 
 *    compute HLL conserved quantities using Mignone eq 9
 * */
   Fl[0] = CL[0] * lV1;
   Fl[1] = CL[1] * lV1 + PL[4];
   Fl[2] = CL[2] * lV1;
   Fl[3] = CL[3] * lV1;
   Fl[4] = CL[1];

   Fr[0] = CR[0] * rV1;
   Fr[1] = CR[1] * rV1 + PR[4];
   Fr[2] = CR[2] * rV1;
   Fr[3] = CR[3] * rV1;
   Fr[4] = CR[1];

 
/* 5. Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 *    Compute HLL conserved quantities using Mignone eq 9
 */
  ovlrmll = 1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;

  Fhll[0] = (lmdar*Fl[0] - lmdal*Fr[0] + lmdatlmda * (CR[0] - CL[0])) * ovlrmll;
  Fhll[1] = (lmdar*Fl[1] - lmdal*Fr[1] + lmdatlmda * (CR[1] - CL[1])) * ovlrmll;
  Fhll[2] = (lmdar*Fl[2] - lmdal*Fr[2] + lmdatlmda * (CR[2] - CL[2])) * ovlrmll;
  Fhll[3] = (lmdar*Fl[3] - lmdal*Fr[3] + lmdatlmda * (CR[3] - CL[3])) * ovlrmll;
  Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] - CL[4])) * ovlrmll;

  Uhll[0] = (lmdar * CR[0] - lmdal * CL[0] + Fl[0] - Fr[0]) * ovlrmll;
  Uhll[1] = (lmdar * CR[1] - lmdal * CL[1] + Fl[1] - Fr[1]) * ovlrmll;
  Uhll[2] = (lmdar * CR[2] - lmdal * CL[2] + Fl[2] - Fr[2]) * ovlrmll;
  Uhll[3] = (lmdar * CR[3] - lmdal * CL[3] + Fl[3] - Fr[3]) * ovlrmll;
  Uhll[4] = (lmdar * CR[4] - lmdal * CL[4] + Fl[4] - Fr[4]) * ovlrmll;

/* 6. Compute contact wave speed using larger root from Mignone Eq 18
 *    Physical root is the root with the minus sign
 */
  /* quadratic formCLa calcCLation */

  a = Fhll[4];
  b = -(Uhll[4] + Fhll[1]);
  c = Uhll[1];

  quad = -0.5*(b + SIGN(b)*SQRT(b*b - 4.0*a*c));
  lmdas = c/quad;

 /* 7. Determine intercell flux according to Mignone 13
 */
  if( lmdal >= 0.0){ /* Fl */
    /* intercell flux is left flux */
    Flux_Out[0] = Fl[0];
    Flux_Out[1] = Fl[1];
    Flux_Out[2] = Fl[2];
    Flux_Out[3] = Fl[3];
    Flux_Out[4] = Fl[4];
   return;
  }
  else if( lmdas >= 0.0){ /* Fls */

    /* Mignone 2006 Eq 48 */
    ps = -Fhll[4]*lmdas + Fhll[1];

    /* now calcCLate Usl with Mignone Eq 16 */
    den = 1.0 / (lmdal - lmdas);

    Usl[0] =  CL[0] * (lmdal - lV1) * den;
    Usl[1] = (CL[1] * (lmdal - lV1) + ps - PL[4]) * den;
    Usl[2] =  CL[2] * (lmdal - lV1) * den;
    Usl[3] =  CL[3] * (lmdal - lV1) * den;
    Usl[4] = (CL[4] * (lmdal - lV1) + ps * lmdas - PL[4] * lV1) * den;

    /* now calcCLate Fsr using Mignone Eq 14 */

    Flux_Out[0] = lmdal*(Usl[0] - CL[0]) + Fl[0];
    Flux_Out[1] = lmdal*(Usl[1] - CL[1]) + Fl[1];
    Flux_Out[2] = lmdal*(Usl[2] - CL[2]) + Fl[2];
    Flux_Out[3] = lmdal*(Usl[3] - CL[3]) + Fl[3];
    Flux_Out[4] = lmdal*(Usl[4] - CL[4]) + Fl[4];

    return;
  }
  else if( lmdar >= 0.0){ /* Frs */
    int count=0;
    if(count==1){
    for(int i=0;i<5;i++)
    printf("PR[%d]=%f, PL[%d]=%f\n",i,PR[i],i,PL[i]);
    abort();
    }
    count++;

    /* Mignone 2006 Eq 48 */
    ps = -Fhll[4]*lmdas + Fhll[1];

    /* now calcCLate Usr with Mignone Eq 16 */
    den = 1.0 / (lmdar - lmdas);

    Usr[0] =  CR[0] * (lmdar - rV1) * den;
    Usr[1] = (CR[1] * (lmdar - rV1) + ps - PR[4]) * den;
    Usr[2] =  CR[2] * (lmdar - rV1) * den;
    Usr[3] =  CR[3] * (lmdar - rV1) * den;
    Usr[4] = (CR[4] * (lmdar - rV1) + ps * lmdas - PR[4] * rV1) * den;

    /* now calcCLate Fsr using Mignone Eq 14 */
    Flux_Out[0] = lmdar*(Usr[0] - CR[0]) + Fr[0];
    Flux_Out[1] = lmdar*(Usr[1] - CR[1]) + Fr[1];
    Flux_Out[2] = lmdar*(Usr[2] - CR[2]) + Fr[2];
    Flux_Out[3] = lmdar*(Usr[3] - CR[3]) + Fr[3];
    Flux_Out[4] = lmdar*(Usr[4] - CR[4]) + Fr[4];

    return;
  }
  else{ /* Fr */
    /* intercell flux is right flux */
    Flux_Out[0] = Fr[0];
    Flux_Out[1] = Fr[1];
    Flux_Out[2] = Fr[2];
    Flux_Out[3] = Fr[3];
    Flux_Out[4] = Fr[4];

    return;
  }

/* 8. restore the correct order */
   CPU_Rotate3D( Flux_Out, XYZ, false );


} // FUNCTION : CPU_RiemannSolver_HLLC



#endif // #if ( MODEL == SR_HYDRO )
