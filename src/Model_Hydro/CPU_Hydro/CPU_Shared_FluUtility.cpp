#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__



#include "CUFLU.h"
#include <stdio.h>

#if ( MODEL == HYDRO )



// internal function prototypes
// --> only necessary for GPU since they are included in Prototype.h for the CPU codes
#ifdef __CUDACC__
GPU_DEVICE
real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinPres, const real MinPres, const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde,
                     const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut );
GPU_DEVICE
real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const bool CheckMinEint, const real MinEint, const real Emag );
GPU_DEVICE
static real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint,
                                const real Emag );
GPU_DEVICE
static real Hydro_CheckMinPres( const real InPres, const real MinPres );

#ifdef SRHD
GPU_DEVICE
real SRHD_Con2HTilde( const real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                      const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                      const real *const EoS_Table[EOS_NTABLE_MAX] );

GPU_DEVICE
void SRHD_HTildeFunction (real HTilde, real MSqr_DSqr, real Temp, real Constant,
                          const EoS_H2TEM_t EoS_HTilde2Temp, real *Fun, real *DiffFun,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] );
GPU_DEVICE
bool SRHD_CheckUnphysical( const real Con[], const real Pri[], const char s[], const int line, bool show );
#endif
#endif

GPU_DEVICE
void  NewtonRaphsonSolver(void (*FunPtr)(real, real, real, real, const EoS_H2TEM_t, real*, real*, const double*, const int*, const real *const*),
                          real MSqr_DSqr, real Constant, real *root, const EoS_H2TEM_t EoS_HTilde2Temp,
                          const real guess, const real epsabs, const real epsrel,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] );
GPU_DEVICE
real VectorDotProduct( real V1, real V2, real V3 );




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Rotate3D
// Description :  Rotate the input fluid variables properly to simplify the 3D calculation
//
// Note        :  1. x : (x,y,z) <--> (x,y,z)
//                   y : (x,y,z) <--> (y,z,x)
//                   z : (x,y,z) <--> (z,x,y)
//                2. Work no matter InOut[] includes passive scalars or not since they are not modified at all
//                   --> For MHD, specify the array offset of magnetic field by Mag_Offset
//
// Parameter   :  InOut      : Array storing both the input and output data
//                XYZ        : Target spatial direction : (0/1/2) --> (x/y/z)
//                Forward    : (true/false) <--> (forward/backward)
//                Mag_Offset : Array offset of magnetic field (for MHD only)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset )
{

   if ( XYZ == 0 )   return;


// check
#  ifdef GAMER_DEBUG
#  ifdef MHD
   if ( Mag_Offset < NCOMP_FLUID  ||  Mag_Offset > NCOMP_TOTAL_PLUS_MAG - NCOMP_MAG )
      printf( "ERROR : invalid Mag_Offset = %d !!\n", Mag_Offset );
#  endif
#  endif


   real Temp_Flu[3];
   for (int v=0; v<3; v++)    Temp_Flu[v] = InOut[ v + 1 ];
#  ifdef MHD
   real Temp_Mag[3];
   for (int v=0; v<3; v++)    Temp_Mag[v] = InOut[ v + Mag_Offset ];
#  endif

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : InOut[              1 ] = Temp_Flu[1];
                  InOut[              2 ] = Temp_Flu[2];
                  InOut[              3 ] = Temp_Flu[0];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[1];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[0];
#                 endif
                  break;

         case 2 : InOut[              1 ] = Temp_Flu[2];
                  InOut[              2 ] = Temp_Flu[0];
                  InOut[              3 ] = Temp_Flu[1];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[0];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[1];
#                 endif
                  break;
      }
   }

   else // backward
   {
      switch ( XYZ )
      {
         case 1 : InOut[              1 ] = Temp_Flu[2];
                  InOut[              2 ] = Temp_Flu[0];
                  InOut[              3 ] = Temp_Flu[1];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[0];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[1];
#                 endif
                  break;

         case 2 : InOut[              1 ] = Temp_Flu[1];
                  InOut[              2 ] = Temp_Flu[2];
                  InOut[              3 ] = Temp_Flu[0];
#                 ifdef MHD
                  InOut[ Mag_Offset + 0 ] = Temp_Mag[1];
                  InOut[ Mag_Offset + 1 ] = Temp_Mag[2];
                  InOut[ Mag_Offset + 2 ] = Temp_Mag[0];
#                 endif
                  break;
      }
   }

} // FUNCTION : Hydro_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Pri
// Description :  Conserved variables --> primitive variables
//
// Note        :  1. Always apply pressure floor
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling Hydro_NormalizePassive() in Hydro_Shared_FullStepUpdate()
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//
// Parameter   :  In                 : Input conserved variables
//                Out                : Output primitive variables
//                MinPres            : Minimum allowed pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS_DensEint2Pres  : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint  : EoS routine to compute the gas internal energy
//                EoS_GuessHTilde    :        . . .               gussed reduced specific enthalpy
//                EoS_HTilde2Temp    :        . . .               temperature
//                EoS_AuxArray_*     : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table          : EoS tables for EoS_DensEint2Pres()
//                EintOut            : Pointer to store the output internal energy
//                                     --> Do nothing if it is NULL
//                                     --> Internal energy floor is not applied
//                LorentzFactor_Ptr  : Pointer to store Lorentz factor
//
// Return      :  Out[], EintOut (optional)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactor_Ptr )
{

#  ifndef SRHD
   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0/In[0];
#  endif

#  if( !defined SRHD && !defined MHD )
   const real Emag             = NULL_REAL;
#  endif

#  ifdef MHD
   const real Bx               = In[ MAG_OFFSET + 0 ];
   const real By               = In[ MAG_OFFSET + 1 ];
   const real Bz               = In[ MAG_OFFSET + 2 ];
   const real Emag             = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else

#  endif

// conserved --> primitive
#  ifdef SRHD
   real HTilde, Factor, Temp, AAA;
   real *LorentzFactor = &AAA;
   SRHD_CheckUnphysical( In, NULL, __FUNCTION__, __LINE__, true  );
   HTilde = SRHD_Con2HTilde( In, EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );


   Factor = In[0]*((real)1.0 + HTilde);

   Out[1] = In[1]/Factor;
   Out[2] = In[2]/Factor;
   Out[3] = In[3]/Factor;

   if ( LorentzFactor_Ptr != NULL )  LorentzFactor = LorentzFactor_Ptr;

   *LorentzFactor = SQRT( (real)1.0 + SQR(Out[1]) + SQR(Out[2]) + SQR(Out[3]) );

   Out[0] = In[0]/ *LorentzFactor;

   EoS_HTilde2Temp( HTilde, &Temp, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Out[4] = Out[0]*Temp;
#  else
   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = Hydro_Con2Pres( In[0], In[1], In[2], In[3], In[4], In+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag,
                            EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, EintOut );


// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   {
      const real Pres0 = Out[4];
      Out[4] = Hydro_CheckMinPres( Pres0, JeansMinPres_Coeff*SQR(Out[0]) );

//    recompute internal energy to be consistent with the updated pressure
      if ( EintOut != NULL  &&  Out[4] != Pres0 )
         *EintOut = EoS_DensPres2Eint( Out[0], Out[4], In+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   }

#  endif

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass density of target passive scalars to mass fraction
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= _Rho;
#  endif


// B field
#  ifdef MHD
   for (int v=NCOMP_TOTAL; v<NCOMP_TOTAL_PLUS_MAG; v++)  Out[v] = In[v];
#  endif

} // FUNCTION : Hydro_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Pri2Con
// Description :  Primitive variables --> conserved variables
//
// Note        :  1. Does NOT check if the input pressure is greater than the given minimum threshold
//                2. For passive scalars, we store their mass fraction as the primitive variables
//                   when NormPassive is on
//                   --> See the input parameters "NormPassive, NNorm, NormIdx"
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//                5. Convert pressure to internal energy using the input EoS routine by default
//                   --> But one can also specify internal energy directly through *EintIn*,
//                       by which no EoS conversion is required and the input pressure will be useless
//                       --> Mainly used by the option LR_EINT in data reconstruction
//
// Parameter   :  In                : Array storing the input primitive variables
//                Out               : Array to store the output conserved variables
//                NormPassive       : true --> convert passive scalars to mass fraction
//                NNorm             : Number of passive scalars for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx           : Target variable indices for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_VarIdx"
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_Temp2HTilde   : EoS routine to compute the reduced specific enthalpy
//                EoS_HTilde2Temp   : EoS routine to compute the temperature
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensPres2Eint()
//                EoS_Table         : EoS tables for EoS_DensPres2Eint()
//                EintIn            : Pointer storing the input internal energy (see the note above)
//                                    --> Do nothing if it is NULL
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Pri2Con( const real In[], real Out[], const bool NormPassive, const int NNorm, const int NormIdx[],
                    const EoS_DP2E_t EoS_DensPres2Eint, const EoS_TEM2H_t EoS_Temp2HTilde,
                    const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn )
{

#  ifdef SRHD
   real LorentzFactor, Temperature, HTilde, MSqr_DSqr, HTildeFunction, Factor0;
#  else
   real Eint;
   real Emag=NULL_REAL;
#  endif


// passive scalars
// --> do it before invoking EoS_DensPres2Eint() since the latter requires the mass density
//     instead of mass fraction of passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( NormPassive )
      for (int v=0; v<NNorm; v++)   Out[ NCOMP_FLUID + NormIdx[v] ] *= In[0];
#  endif


// primitive --> conserved
#  ifdef SRHD
   LorentzFactor = SQRT( (real)1.0 + SQR(In[1]) + SQR(In[2]) + SQR(In[3]) );
   Temperature = In[4]/In[0];
   HTilde = EoS_Temp2HTilde( Temperature, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Out[0] = In[0]*LorentzFactor;
   Factor0 = Out[0]*HTilde + Out[0];

   Out[1] = In[1]*Factor0;
   Out[2] = In[2]*Factor0;
   Out[3] = In[3]*Factor0;
   MSqr_DSqr  = SQR(Out[1])+SQR(Out[2])+SQR(Out[3]);
   MSqr_DSqr /= SQR(Out[0]);
   SRHD_HTildeFunction( HTilde, MSqr_DSqr, Temperature, (real)0.0, EoS_HTilde2Temp, &HTildeFunction, NULL,
                        EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Out[4]  = MSqr_DSqr + HTildeFunction;
   Out[4] /= (real)1.0 + SQRT( (real)1.0 + MSqr_DSqr + HTildeFunction );
   Out[4] *= Out[0];
   SRHD_CheckUnphysical( Out, NULL, __FUNCTION__, __LINE__, true  );
#  else
   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
#  endif

#  ifdef MHD
   const real Bx = In[ MAG_OFFSET + 0 ];
   const real By = In[ MAG_OFFSET + 1 ];
   const real Bz = In[ MAG_OFFSET + 2 ];
   Emag   = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  endif

#  ifndef SRHD
   Eint   = ( EintIn == NULL ) ? EoS_DensPres2Eint( In[0], In[4], Out+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                    EoS_AuxArray_Int, EoS_Table )
                               : *EintIn;
   Out[4] = Hydro_ConEint2Etot( Out[0], Out[1], Out[2], Out[3], Eint, Emag );
#  endif


// B field
#  ifdef MHD
   for (int v=NCOMP_TOTAL; v<NCOMP_TOTAL_PLUS_MAG; v++)  Out[v] = In[v];
#  endif

} // FUNCTION : Hydro_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Flux
// Description :  Evaluate hydrodynamic/MHD fluxes from the input conserved variables
//
// Note        :  1. Flux[] and In[] may point to the same array
//                2. Flux[] and In[] should have the size of NCOMP_TOTAL_PLUS_MAG
//                3. By default, it computes pressure using the input EoS routine
//                   --> But one can also specify pressure directly through *PresIn*,
//                       by which no EoS conversion is required
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux              : Array to store the output fluxes
//                In                : Array storing the input conserved variables
//                MinPres           : Minimum allowed pressure
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables for EoS_DensEint2Pres()
//                AuxArray          : Pointer storing
//                                   (1) the input pressure in non-SRHD (see the note above)
//                                    --> Do nothing if it is NULL
//                                   (2) primitive variables in SRHD
//
// Return      :  Flux[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real AuxArray[] )
{

#  ifndef SRHD
   const bool CheckMinPres_Yes = true;
#  endif

   real InRot[ NCOMP_FLUID + NCOMP_MAG ];    // no need to include passive scalars since they don't have to be rotated
#  ifdef SRHD
   real PriRot[ NCOMP_FLUID + NCOMP_MAG ];
#  endif

   for (int v=0; v<NCOMP_FLUID; v++)   InRot[v] = In[v];
#  ifdef SRHD
   for (int v=0; v<NCOMP_FLUID; v++)   PriRot[v] = AuxArray[v];
#  endif

#  ifdef MHD
   for (int v=NCOMP_FLUID; v<NCOMP_FLUID+NCOMP_MAG; v++)    InRot[v] = In[ v - NCOMP_FLUID + MAG_OFFSET ];
#  endif

   Hydro_Rotate3D( InRot, XYZ, true, NCOMP_FLUID );
#  ifdef SRHD
   Hydro_Rotate3D( PriRot, XYZ, true, NCOMP_FLUID );
#  endif

#  ifdef MHD
   const real Bx   = InRot[ NCOMP_FLUID + 0 ];
   const real By   = InRot[ NCOMP_FLUID + 1 ];
   const real Bz   = InRot[ NCOMP_FLUID + 2 ];
   const real Emag = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  elif ( !defined SRHD )
   const real Emag = NULL_REAL;
#  endif

#  ifdef SRHD
   const real LorentzFactor = SQRT((real)1.0 + VectorDotProduct(PriRot[1], PriRot[2], PriRot[3]));
   const real Vx = PriRot[1] / LorentzFactor;

   Flux[0] = InRot[0]*Vx;
   Flux[1] = FMA( InRot[1], Vx, PriRot[4] );
   Flux[2] = InRot[2] * Vx;
   Flux[3] = InRot[3] * Vx;
   Flux[4] = ( InRot[4] + PriRot[4] )*Vx;
#  else
   const real Pres = ( AuxArray == NULL ) ? Hydro_Con2Pres( InRot[0], InRot[1], InRot[2], InRot[3], InRot[4], In+NCOMP_FLUID,
                                                            CheckMinPres_Yes, MinPres, Emag,
                                                            EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                            EoS_Table, NULL )
                                        : AuxArray[0];
   const real _Rho = (real)1.0 / InRot[0];
   const real Vx   = _Rho*InRot[1];

   Flux[0] = InRot[1];
   Flux[1] = Vx*InRot[1] + Pres;
   Flux[2] = Vx*InRot[2];
   Flux[3] = Vx*InRot[3];
   Flux[4] = Vx*( InRot[4] + Pres );
#  endif

// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux[v] = In[v]*Vx;
#  endif

// B field
#  ifdef MHD
   const real Vy = _Rho*InRot[2];
   const real Vz = _Rho*InRot[3];

   Flux[              1 ] += Emag - SQR(Bx);
   Flux[              2 ] -= Bx*By;
   Flux[              3 ] -= Bx*Bz;
   Flux[              4 ] += Vx*Emag - Bx*( Bx*Vx + By*Vy + Bz*Vz );
   Flux[ MAG_OFFSET + 0 ]  = (real)0.0;
   Flux[ MAG_OFFSET + 1 ]  = By*Vx - Bx*Vy;
   Flux[ MAG_OFFSET + 2 ]  = Bz*Vx - Bx*Vz;
#  endif

   Hydro_Rotate3D( Flux, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_Con2Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHD_Con2HTilde
// Description :  
// Note        :  
// Parameter   :  Con             : Array storing the conservative variables
//                EoS_GuessHTilde : EoS routine to compute the reduced energy
//                EoS_HTilde2Temp : EoS routine to compute the temperature
//-------------------------------------------------------------------------------------------------------
#ifdef SRHD
GPU_DEVICE
real SRHD_Con2HTilde( const real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                      const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                      const real *const EoS_Table[EOS_NTABLE_MAX] )
{
  real HTilde, GuessHTilde, MSqr_DSqr, Constant;

  MSqr_DSqr  = SQR(Con[1])+SQR(Con[2])+SQR(Con[3]);
  MSqr_DSqr /= SQR(Con[0]);

  GuessHTilde = EoS_GuessHTilde( Con, &Constant, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table);

  void (*FunPtr)( real HTilde, real MSqr_DSqr, real Temp, real Constant,
                  const EoS_H2TEM_t EoS_HTilde2Temp, real *Fun, real *DiffFun,
                  const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                  const real *const EoS_Table[EOS_NTABLE_MAX] ) = &SRHD_HTildeFunction;

  NewtonRaphsonSolver(FunPtr, MSqr_DSqr, -Constant, &HTilde, EoS_HTilde2Temp, GuessHTilde,
                      (real)TINY_NUMBER, (real)MACHINE_EPSILON, EoS_AuxArray_Flt,
                      EoS_AuxArray_Int, EoS_Table );

  return HTilde;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHD_HTildeFunction
// Description :  
// Note        :  
// Parameter   :  HTilde          : The reduced specific enthalpy
//                MSqr_DSqr       : (|Momentum|/Dens)**2
//                Temp            : The temperature
//                Constant        : The constant on the other side of the function to be iterated
//                EoS_HTilde2Temp : EoS routine to compute the temperature
//               *Fun             : The function to be numerically solved
//               *DiffFun         : The derivative function with respect to the unknown variable
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHD_HTildeFunction (real HTilde, real MSqr_DSqr, real Temp, real Constant,
                          const EoS_H2TEM_t EoS_HTilde2Temp, real *Fun, real *DiffFun,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{
  real DiffTemp;

  if ( Temp != Temp )
  EoS_HTilde2Temp( HTilde, &Temp, &DiffTemp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );


  real H =  HTilde + (real)1.0;
  real Factor0 = SQR( H ) + MSqr_DSqr;

  if ( Fun != NULL )

  *Fun = SQR( HTilde ) + (real)2.0*HTilde - (real)2.0*Temp - (real)2.0*Temp*HTilde
		  + SQR( Temp * H ) / Factor0 + Constant;

  if ( DiffFun != NULL )

  *DiffFun = (real)2.0*H - (real)2.0*Temp - (real)2.0*H*DiffTemp +
		  ( (real)2.0*Temp*DiffTemp*H*H - (real)2.0*Temp*Temp*H ) / SQR( Factor0 );

}

GPU_DEVICE
real SRHD_Con2KineticEngy( real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                           const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                           const real *const EoS_Table[EOS_NTABLE_MAX] )
{
  real H, Usqr, Pri[NCOMP_FLUID], LorentzFactor;

  H = (real)1.0 + SRHD_Con2HTilde( Con, EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

  Hydro_Con2Pri( Con, Pri, NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL, NULL, NULL,
                 EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, &LorentzFactor );

  Usqr = VectorDotProduct( Pri[1], Pri[2], Pri[3]  );

  return ( Con[DENS] * H + Pri[4] ) * Usqr / ( LorentzFactor + (real)1.0 );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHD_CheckUnphysical
// Description :  
// Note        :  
// Parameter   :  Con           : Array storing the conservative variables
//                Pri           : Array storing the primitive variables
//                FunctionName  : The function occurring unphysical result
//                Line          : The line in `FunctionName`
//                Show          : Print unphysical result in log or not
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool SRHD_CheckUnphysical( const real Con[], const real Pri[], const char FunctionName[], const int Line, bool Show )
{
   real discriminant;
   real Msqr, M, E_D, M_D;


//--------------------------------------------------------------//
//------------ only check conserved variables-------------------//
//--------------------------------------------------------------//
   if ( Con != NULL && Pri == NULL)
   {
// check NaN
      if (  Con[DENS] != Con[DENS]
         || Con[MOMX] != Con[MOMX]
         || Con[MOMY] != Con[MOMY]
         || Con[MOMZ] != Con[MOMZ]
         || Con[ENGY] != Con[ENGY]  )                                                   goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= Con[DENS] || Con[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMX] || Con[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMY] || Con[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Con[MOMZ] || Con[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Con[ENGY] || Con[ENGY]  >= (real)HUGE_NUMBER )       goto FAIL;


// check minimum energy
      Msqr = VectorDotProduct( Con[MOMX], Con[MOMY], Con[MOMZ] );
      M = SQRT( Msqr );
	  E_D = Con[ENGY] / Con[DENS];
	  M_D = M / Con[DENS];
#     ifdef REDUCED_ENERGY
      discriminant = ( E_D + M_D ) * ( E_D - M_D ) + (real)2.0 * E_D;
      if ( discriminant <= TINY_NUMBER )                                                goto FAIL;
#     else
      discriminant = ( ( SQR( Con[ENGY] ) -  Msqr ) / SQR ( Con[DENS] ) );
      if ( discriminant <= 1.0   )                                                      goto FAIL;
#     endif



// pass all checks
      return false;
   }

//--------------------------------------------------------------//
//------------ only check primitive variables-------------------//
//--------------------------------------------------------------//

   else if ( Con == NULL && Pri != NULL)
   {
// check NaN
      if (  Pri[DENS] != Pri[DENS]
         || Pri[MOMX] != Pri[MOMX]
         || Pri[MOMY] != Pri[MOMY]
         || Pri[MOMZ] != Pri[MOMZ]
         || Pri[ENGY] != Pri[ENGY]  )                                                   goto FAIL;

// check +inf and -inf
      if (  (real)  TINY_NUMBER >= Pri[DENS] || Pri[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMX] || Pri[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMY] || Pri[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Pri[MOMZ] || Pri[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Pri[ENGY] || Pri[ENGY]  >= (real)HUGE_NUMBER )       goto FAIL;


// pass all checks
      return false;
   }

// print all variables if goto FAIL
      FAIL:
      {
        if ( Show )
         {
           printf( "\n\nError!! function: %s: %d\n", FunctionName, Line);

           if ( Con != NULL && Pri == NULL)
           {
              printf( "D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n",
                                   Con[DENS], Con[MOMX], Con[MOMY], Con[MOMZ], Con[ENGY]);
              printf( "E^2+2*E*D-|M|^2=%14.7e\n", discriminant );
           }
           else
           {
              printf( "n=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n",
                                   Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);
           }

         }
        return true;
      }
}
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinPres
// Description :  Check if the input pressure is great than the minimum allowed threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high Mach numbers
//                   --> Currently it simply sets a minimum allowed value for pressure
//                       --> Set MIN_PRES in the runtime parameter file "Input__Parameter"
//                2. We should also support a minimum **temperature** instead of **pressure**
//                   --> NOT supported yet
//                3. If the input pressure is NaN, return NaN in order to trigger auto-correction such as
//                   "OPT__1ST_FLUX_CORR" and "AUTO_REDUCE_DT"
//
// Parameter   :  InPres  : Input pressure to be corrected
//                MinPres : Minimum allowed pressure
//
// Return      :  InPres != NaN --> max( InPres, MinPres )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinPres( const real InPres, const real MinPres )
{

// call FMAX() only if InPres is not NaN
   if ( InPres == InPres )    return FMAX( InPres, MinPres );
   else                       return InPres;

} // FUNCTION : Hydro_CheckMinPres



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinEint
// Description :  Similar to Hydro_CheckMinPres() except that this function checks the internal energy
//                density (Eint) instead of pressure
//
// Note        :  1. See Hydro_CheckMinPres()
//
// Parameter   :  InEint  : Input Eint to be corrected
//                MinEint : Minimum allowed Eint
//
// Return      :  InEint != NaN --> max( InEint, MinEint )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEint( const real InEint, const real MinEint )
{

// call FMAX() only if InEint is not NaN
   if ( InEint == InEint )    return FMAX( InEint, MinEint );
   else                       return InEint;

} // FUNCTION : Hydro_CheckMinEint



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckMinEintInEngy
// Description :  Ensure that the internal energy density in the input total energy density is greater than
//                a given threshold
//
// Note        :  1. Invoke Hydro_CheckMinEint()
//                2. Input conserved instead of primitive variables
//                3. For MHD, one must provide the magnetic energy density Emag (i.e., 0.5*B^2)
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                InEngy   : Energy density
//                MinEint  : Internal energy density floor
//                Emag    : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Total energy density with internal energy density greater than a given threshold
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEintInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real InEngy,
                               const real MinEint, const real Emag )
{

   const bool CheckMinEint_No = false;
   real InEint, OutEint, OutEngy;

   InEint  = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, InEngy, NULL, NULL, CheckMinEint_No, NULL_REAL, Emag );
   OutEint = Hydro_CheckMinEint( InEint, MinEint );

// do not modify energy (even the round-off errors) if the input data pass the check
   if ( InEint == OutEint )   OutEngy = InEngy;
   else                       OutEngy = InEngy - InEint + OutEint;

   return OutEngy;

} // FUNCTION : Hydro_CheckMinEintInEngy



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckNegative
// Description :  Check whether the input value is <= 0.0 (also check whether it's Inf or NAN)
//
// Note        :  Can be used to check whether the values of density and pressure are unphysical
//
// Parameter   :  Input : Input value
//
// Return      :  true  --> Input <= 0.0  ||  >= __FLT_MAX__  ||  != itself (Nan)
//                false --> otherwise
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool Hydro_CheckNegative( const real Input )
{

   if ( Input <= (real)0.0  ||  Input >= __FLT_MAX__  ||  Input != Input )    return true;
   else                                                                       return false;

} // FUNCTION : Hydro_CheckNegative



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_ConEint2Etot
// Description :  Evaluate total energy from the input conserved variables and internal energy
//
// Note        :  1. For MHD, total energy density includes the magnetic energy Emag=0.5*B^2
//                2. Internal energy density is energy per volume instead of per mass
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Eint     : Internal energy density
//                Emag     : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Total energy density (including the magnetic energy density for MHD)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint,
                         const real Emag )
{

//###NOTE: assuming Etot = Eint + Ekin + Emag
   real Etot;

   Etot  = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / Dens;
   Etot += Eint;
#  ifdef MHD
   Etot += Emag;
#  endif

   return Etot;

} // FUNCTION : Hydro_ConEint2Etot



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Temp2Pres
// Description :  Convert gas temperature to pressure
//
// Note        :  1. Assume the ideal-gas law
//                   --> P = \rho*K*T / ( mu*m_H )
//                2. Assume both input and output to be code units
//                   --> Temperature should be converted to UNIT_E in advance
//                       --> Example: T_code_unit = T_kelvin * Const_kB / UNIT_E
//                3. Pressure floor (MinPres) is applied when enabling CheckMinPres
//                4. Adopt double precision since
//                   (1) both Temp and m_H may exhibit extreme values depending on the code units, and
//                   (2) we don't really care about the performance here since this function is usually
//                       only used for constructing the initial condition
//
// Parameter   :  Dens         : Gas mass density in code units
//                Temp         : Gas temperature in code units
//                mu           : Mean molecular weight
//                m_H          : Atomic hydrogen mass in code units
//                               --> Sometimes we use the atomic mass unit (Const_amu defined in PhysicalConstant.h)
//                                   and m_H (Const_mH defined in PhysicalConstant.h) interchangeably since the
//                                   difference is small (m_H ~ 1.007825 amu)
//                CheckMinPres : Apply pressure floor by calling Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Pressure floor
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
double Hydro_Temp2Pres( const double Dens, const double Temp, const double mu, const double m_H,
                        const bool CheckMinPres, const double MinPres )
{

   double Pres;

   Pres = Dens*Temp/(mu*m_H);

   if ( CheckMinPres )  Pres = Hydro_CheckMinPres( (real)Pres, (real)MinPres );

   return Pres;

} // FUNCTION : Hydro_Temp2Pres




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Pres
// Description :  Evaluate the fluid pressure
//
// Note        :  1. Invoke the EoS routine EoS_DensEint2Pres() to support different EoS
//                2. For MHD, Engy is the total energy density including the magnetic energy Emag=0.5*B^2
//                   and thus one must provide Emag to subtract it
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density (including the magnetic energy density for MHD)
//                Passive           : Passive scalars
//                CheckMinPres      : Apply pressure floor by invoking Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which this option should be disabled
//                                        --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres           : Pressure floor
//                Emag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables for EoS_DensEint2Pres()
//                EintOut           : Pointer to store the output internal energy
//                                    --> Do nothing if it is NULL
//                                    --> Internal energy floor is not applied
//
// Return      :  Gas pressure (Pres), EintOut (optional)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinPres, const real MinPres, const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde,
                     const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut )
{

   real Pres;

#  ifdef SRHD
   real Cons[NCOMP_FLUID] = { Dens, MomX, MomY, MomZ, Engy };
   real Prim[NCOMP_FLUID];
   Hydro_Con2Pri( Cons, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, (real)NULL_REAL,
                  NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, NULL );
   Pres = Prim[4];
#  else
   const bool CheckMinEint_No = false;
   real Eint;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, NULL, NULL, CheckMinEint_No, NULL_REAL, Emag );
   Pres = EoS_DensEint2Pres( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

   if ( CheckMinPres )   Pres = Hydro_CheckMinPres( Pres, MinPres );

   if ( EintOut != NULL )  *EintOut = Eint;
#  endif

   return Pres;

} // FUNCTION : Hydro_Con2Pres



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Eint
// Description :  Evaluate the gas internal energy density
//
// Note        :  1. For MHD, Engy is the total energy density including the magnetic energy Emag=0.5*B^2
//                   and thus one must provide Emag to subtract it
//                2. Internal energy density is energy per volume instead of per mass
//
// Parameter   :  Dens         : Mass density
//                MomX/Y/Z     : Momentum density
//                Engy         : Energy density (including the magnetic energy density for MHD)
//                CheckMinEint : Apply internal energy floor by invoking Hydro_CheckMinEint()
//                               --> In some cases we actually want to check if internal energy becomes unphysical,
//                                   for which this option should be disabled
//                MinEint      : Internal energy floor
//                Emag         : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Gas internal energy density (Eint)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const bool CheckMinEint, const real MinEint, const real Emag )
{

//###NOTE: assuming Etot = Eint + Ekin + Emag
   real Eint;

#  ifdef SRHD
   real Cons[NCOMP_FLUID] = { Dens, MomX, MomY, MomZ, Engy };
   real Prim[NCOMP_FLUID], HTilde_1, Temp, TempSqr;

   printf( "%s: Hydro_Con2Pri need EoS table !!\n", __FUNCTION__);
   exit(0);

   Hydro_Con2Pri( Cons, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                  (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  NULL, NULL, NULL, NULL, NULL );

   Temp = Prim[4]/Prim[0];
   TempSqr =  SQR(Temp);
   HTilde_1  = (real)2.25 * TempSqr;
   HTilde_1 /= (real)1.0 + SQRT( (real)1.0 + (real)2.25*TempSqr );
   HTilde_1 += (real)1.5 * Temp;
   Eint    = Prim[0] * HTilde_1;
#  else
   Eint  = Engy - (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / Dens;
#  endif

#  ifdef MHD
   Eint -= Emag;
#  endif

#  if ( !defined SRHD && defined MHD )
   if ( CheckMinEint )   Eint = Hydro_CheckMinEint( Eint, MinEint );
#  endif

   return Eint;

} // FUNCTION : Hydro_Con2Eint



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Temp
// Description :  Evaluate the fluid temperature
//
// Note        :  1. For simplicity, currently this function only returns **pressure/density**, which does
//                   NOT include normalization
//                   --> For OPT__FLAG_LOHNER_TEMP only
//                   --> Also note that currently it only checks minimum pressure but not minimum density
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density
//                Passive           : Passive scalars
//                CheckMinPres      : Apply pressure floor by calling Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which we don't want to enable this option
//                                        --> For example: Flu_FixUp(), Flu_Close(), Hydro_Aux_Check_Negative()
//                MinPres           : Pressure floor
//                Bmag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table         : EoS tables for EoS_DensEint2Pres()
//
// Return      :  Gas temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Temp( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinPres, const real MinPres, const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde,
                     const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real Temp;

#  ifdef SRHD
   real Cons[NCOMP_FLUID] =  {Dens, MomX, MomY, MomZ, Engy};
   real Prim[NCOMP_FLUID];
   Hydro_Con2Pri( Cons, Prim, (real)NULL_REAL, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                  (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde,
                  EoS_HTilde2Temp, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, NULL );
   Temp = Prim[4]/Prim[0];
#  else
   const real Pres = Hydro_Con2Pres( Dens, MomX, MomY, MomZ, Engy, Passive, CheckMinPres, MinPres, Emag,
                                     EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

   Temp = Pres/Dens;
#  endif

   return Temp;

} // FUNCTION : Hydro_Con2Temp






//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_NormalizePassive
// Description :  Normalize the target passive scalars so that the sum of their mass density is equal to
//                the gas mass density
//
// Note        :  1. Should be invoked AFTER applying the floor values to passive scalars
//                2. Invoked by Hydro_Shared_FullStepUpdate(), Prepare_PatchData(), Refine(), LB_Refine_AllocateNewPatch(),
//                   Flu_FixUp(), XXX_Init_ByFunction_AssignData(), Flu_Close()
//
// Parameter   :  GasDens : Gas mass density
//                Passive : Passive scalar array (with the size NCOMP_PASSIVE)
//                NNorm   : Number of passive scalars to be normalized
//                          --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx : Target variable indices to be normalized
//                          --> Should be set to the global variable "PassiveNorm_VarIdx"
//
// Return      :  Passive
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] )
{

// validate the target variable indices
#  ifdef GAMER_DEBUG
   const int MinIdx = 0;
   const int MaxIdx = NCOMP_PASSIVE - 1;

   for (int v=0; v<NNorm; v++)
   {
      if ( NormIdx[v] < MinIdx  ||  NormIdx[v] > MaxIdx )
         printf( "ERROR : NormIdx[%d] = %d is not within the correct range ([%d <= idx <= %d]) !!\n",
                 v, NormIdx[v], MinIdx, MaxIdx );
   }
#  endif


   real Norm, PassiveDens_Sum=(real)0.0;

   for (int v=0; v<NNorm; v++)   PassiveDens_Sum += Passive[ NormIdx[v] ];

   Norm = GasDens / PassiveDens_Sum;

   for (int v=0; v<NNorm; v++)   Passive[ NormIdx[v] ] *= Norm;

} // FUNCTION : Hydro_NormalizePassive




#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBField
// Description :  Calculate the cell-centered magnetic field from the input face-centered magnetic field array
//
// Note        :  1. Use the central average operator
//                2. Return all three components of the B field
//                3. Input arrays should have the following dimension:
//                      Bx_FC[]: (Nx+1)*(Ny  )*(Nz  )
//                      By_FC[]: (Nx  )*(Ny+1)*(Nz  )
//                      Bz_FC[]: (Nx  )*(Ny  )*(Nz+1)
//
// Parameter   :  B_CC      : Cell-centered B field to be returned
//                Bx/y/z_FC : Input face-centered B field array
//                Nx/y/z    : Array dimension along different directions (see Note above)
//                i/j/k     : Target cell indices
//
// Return      :  B_CC
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_GetCellCenteredBField( real B_CC[], const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                const int Nx, const int Ny, const int Nz, const int i, const int j, const int k )
{

   const int idx_Bx = IDX321_BX( i, j, k, Nx, Ny );
   const int idx_By = IDX321_BY( i, j, k, Nx, Ny );
   const int idx_Bz = IDX321_BZ( i, j, k, Nx, Ny );

   B_CC[0] = (real)0.5*( Bx_FC[idx_Bx] + Bx_FC[ idx_Bx + 1     ] );
   B_CC[1] = (real)0.5*( By_FC[idx_By] + By_FC[ idx_By + Nx    ] );
   B_CC[2] = (real)0.5*( Bz_FC[idx_Bz] + Bz_FC[ idx_Bz + Nx*Ny ] );

} // FUNCTION : MHD_GetCellCenteredBField



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_GetCellCenteredBEnergy
// Description :  Calculate the cell-centered magnetic energy (i.e., 0.5*B^2) from the input face-centered
//                magnetic field array
//
// Note        :  1. Invoke MHD_GetCellCenteredBField()
//                2. Input arrays should have the following dimension:
//                      Bx_FC[]: (Nx+1)*(Ny  )*(Nz  )
//                      By_FC[]: (Nx  )*(Ny+1)*(Nz  )
//                      Bz_FC[]: (Nx  )*(Ny  )*(Nz+1)
//
// Parameter   :  Bx/y/z_FC : Input face-centered B field array
//                Nx/y/z    : Array dimension along different directions (see Note above)
//                i/j/k     : Target cell indices
//
// Return      :  0.5*B^2 at the center of the target cell
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real MHD_GetCellCenteredBEnergy( const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                 const int Nx, const int Ny, const int Nz, const int i, const int j, const int k )
{

// CC = cell-centered
   real B_CC[3], BEngy;

   MHD_GetCellCenteredBField( B_CC, Bx_FC, By_FC, Bz_FC, Nx, Ny, Nz, i, j, k );

   BEngy = (real)0.5*( SQR(B_CC[MAGX]) + SQR(B_CC[MAGY]) + SQR(B_CC[MAGZ]) );

   return BEngy;

} // FUNCTION : MHD_GetCellCenteredBEnergy
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
// Note        :  
// Parameter   :  
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void  NewtonRaphsonSolver(void (*FunPtr)(real, real, real, real, const EoS_H2TEM_t, real*, real*, const double*, const int*, const real *const*), 
                          real MSqr_DSqr, real Constant, real *root, const EoS_H2TEM_t EoS_HTilde2Temp, 
                          const real guess, const real epsabs, const real epsrel,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{
 int iter = 0;

 #ifdef FLOAT8
 int max_iter = 20;
 #else
 int max_iter = 10;
 #endif

 real Fun, DiffFun;
 real delta;
 real tolerance;
 *root = guess;

 do
   {
     iter++;
     FunPtr(*root, MSqr_DSqr, NAN, Constant, EoS_HTilde2Temp, &Fun, &DiffFun, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

#    ifdef CHECK_FAILED_CELL_IN_FLUID
     if ( DiffFun == (real)0.0 )                                                  printf("derivative is zero\n");
     if ( Fun != Fun  ||(real) -HUGE_NUMBER >= Fun  || Fun  >= (real)HUGE_NUMBER )  printf("function value is not finite\n");
     if ( DiffFun != DiffFun ||(real) -HUGE_NUMBER >= DiffFun || DiffFun >= (real)HUGE_NUMBER )  printf("derivative value is not finite\n");
#    endif

      delta = Fun/DiffFun;
      *root = *root - delta;

      tolerance =  FMA( epsrel, FABS(*root), epsabs );

   }while ( fabs(delta) >= tolerance && iter < max_iter );

}



//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
// Note        :  
// Parameter   :  
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real VectorDotProduct( real V1, real V2, real V3 )
{
  real Product = (real)0.0;

  Product = FMA( V1, V1 , Product );
  Product = FMA( V2, V2 , Product );
  Product = FMA( V3, V3 , Product );

  return Product;
}

#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_FLUUTILITY__
