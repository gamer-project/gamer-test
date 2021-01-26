#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"

void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactor_Ptr );


void Hydro_Pri2Con( const real In[], real Out[], const bool NormPassive, const int NNorm, const int NormIdx[],
                    const EoS_DP2E_t EoS_DensPres2Eint, const EoS_TEM2H_t EoS_Temp2HTilde,
                    const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn );

#ifdef FLOAT8
typedef double real;
#else
typedef float real;
#endif


EoS_GUESS_t  EoS_GuessHTilde_CPUPtr = NULL; 
EoS_H2TEM_t  EoS_HTilde2Temp_CPUPtr = NULL;
EoS_TEM2H_t  EoS_Temp2HTilde_CPUPtr = NULL;
EoS_TEM2C_t  EoS_Temper2CSqr_CPUPtr = NULL;
EoS_DE2P_t   EoS_DensEint2Pres_CPUPtr = NULL; 
EoS_DE2P_t   EoS_DensPres2Eint_CPUPtr = NULL;

double EoS_AuxArray_Flt[EOS_NAUX_MAX];
int EoS_AuxArray_Int[EOS_NAUX_MAX];
real* h_EoS_Table[EOS_NAUX_MAX] = {0};

int MPI_Rank = 0;

int
main ()
{
  EoS_Init();

  real LorentzFactor;

  real Pri[5] = { 0 };
  real Con[5] = { 0 };

  Con[0] = 1.0;
  Con[1] = 0.0;
  Con[2] = 0.0;
  Con[3] = 0.0;
  Con[4] = 1.0;


  printf ("\ninitial conservative vars:\n");
  printf ("=========================================\n\n");

  printf ("D     = %30.25e\n"  ,  Con[0]);
  printf ("Mx    = %30.25e\n"  ,  Con[1]);
  printf ("My    = %30.25e\n"  ,  Con[2]);
  printf ("Mz    = %30.25e\n"  ,  Con[3]);
  printf ("E     = %30.25e\n"  ,  Con[4]);


  real M_Dsqr = SQR( Con[1]/Con[0] ) + SQR( Con[2]/Con[0] ) + SQR( Con[3]/Con[0] );
  real M_D    = SQRT( M_Dsqr );
  real E_D    = Con[4]/Con[0];

  real Discriminant = ( E_D + M_D )*( E_D - M_D ) + (real)2.0*E_D;

  if ( Discriminant < (real)TINY_NUMBER )
  {
		  printf("\nERROR!! Discriminant = %30.25e\n", Discriminant);
		  exit(0);
  }
  else
  {
		  printf("\nDiscriminant = %30.25e\n", Discriminant);
  }

  Hydro_Con2Pri( Con, Pri, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                 (real)NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                 EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL, &LorentzFactor );

  printf ("\nconservative --> primitive\n");
  printf ("=========================================\n\n");

  printf ("n     = %30.25e\n"  ,  Pri[0]);
  printf ("Ux    = %30.25e\n"  ,  Pri[1]);
  printf ("Uy    = %30.25e\n"  ,  Pri[2]);
  printf ("Uz    = %30.25e\n"  ,  Pri[3]);
  printf ("P     = %30.25e\n"  ,  Pri[4]);
  printf ("T     = %30.25e\n"  ,  Pri[4]/Pri[0]);
  printf ("gamma = %30.25e\n"  ,  LorentzFactor);



  real Con_re[5] = { 0 };

  Hydro_Pri2Con( Pri, Con_re, NULL_BOOL, NULL_INT, NULL, NULL,
                 EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );


  printf ("\nprimitive --> conservative\n");
  printf ("=========================================\n\n");
  
  printf ("D     = %30.25e\n"  ,  Con_re[0]);
  printf ("Mx    = %30.25e\n"  ,  Con_re[1]);
  printf ("My    = %30.25e\n"  ,  Con_re[2]);
  printf ("Mz    = %30.25e\n"  ,  Con_re[3]);
  printf ("E     = %30.25e\n"  ,  Con_re[4]);

  real err_D  = (Con_re[0] - Con[0]) / Con[0];
  real err_M1 = (Con_re[1] - Con[1]) / Con[1];
  real err_M2 = (Con_re[2] - Con[2]) / Con[2];
  real err_M3 = (Con_re[3] - Con[3]) / Con[3];
  real err_E  = (Con_re[4] - Con[4]) / Con[4];

  printf ("\nrelative error:\n");
  printf ("=========================================\n\n");

  printf ("err_D   = %+30.25e\n", err_D );

  if ( fabs (Con_re[1]) > TINY_NUMBER ) printf ("err_Mx  = %+30.25e\n", err_M1) ;
  if ( fabs (Con_re[2]) > TINY_NUMBER ) printf ("err_My  = %+30.25e\n", err_M2) ;
  if ( fabs (Con_re[3]) > TINY_NUMBER ) printf ("err_Mz  = %+30.25e\n", err_M3) ;

  printf ("err_E   = %+30.25e\n", err_E );


  return 0;
}
