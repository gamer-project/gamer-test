#ifndef __CUFLU_FLUUTILITY__
#define __CUFLU_FLUUTILITY__



#include "CUFLU.h"
#include "stdio.h"

#if ( MODEL == HYDRO )



// internal function prototypes
// --> only necessary for GPU since they are included in Prototype.h for the CPU codes
#ifdef __CUDACC__
GPU_DEVICE
static real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                            const real Passive[], const bool CheckMinPres, const real MinPres, const real Emag,
                            const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                            const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut );
GPU_DEVICE
static real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                            const bool CheckMinEint, const real MinEint, const real Emag );
GPU_DEVICE
static real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint,
                                const real Emag );
GPU_DEVICE
static real Hydro_CheckMinPres( const real InPres, const real MinPres );
GPU_DEVICE
static real Hydro_CheckMinEint( const real InEint, const real MinEint );
GPU_DEVICE
static real Hydro_CheckMinTemp( const real InTemp, const real MinTemp );
#endif




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
//                   when FracPassive is on
//                   --> See the input parameters "FracPassive, NFrac, FracIdx"
//                   --> But note that here we do NOT ensure "sum(mass fraction) == 1.0"
//                       --> It is done by calling Hydro_NormalizePassive() in Hydro_Shared_FullStepUpdate()
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//
// Parameter   :  In                 : Input conserved variables
//                Out                : Output primitive variables
//                MinPres            : Minimum allowed pressure
//                FracPassive        : true --> convert passive scalars to mass fraction
//                NFrac              : Number of passive scalars for the option "FracPassive"
//                FracIdx            : Target variable indices for the option "FracPassive"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS_DensEint2Pres  : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint  : EoS routine to compute the gas internal energy
//                EoS_AuxArray_*     : Auxiliary arrays for EoS_DensEint2Pres()
//                EoS_Table          : EoS tables for EoS_DensEint2Pres()
//                EintOut            : Pointer to store the output internal energy
//                                     --> Do nothing if it is NULL
//                                     --> Internal energy floor is not applied
//
// Return      :  Out[], EintOut (optional)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut )
{

   const bool CheckMinPres_Yes = true;
   const real _Rho             = (real)1.0/In[0];
#  ifdef MHD
   const real Bx               = In[ MAG_OFFSET + 0 ];
   const real By               = In[ MAG_OFFSET + 1 ];
   const real Bz               = In[ MAG_OFFSET + 2 ];
   const real Emag             = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else
   const real Emag             = NULL_REAL;
#  endif

// conserved --> primitive
   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = Hydro_Con2Pres( In[0], In[1], In[2], In[3], In[4], In+NCOMP_FLUID, CheckMinPres_Yes, MinPres, Emag,
                            EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, EintOut );


// pressure floor required to resolve the Jeans length
// --> note that currently we do not modify the dual-energy variable (e.g., entropy) accordingly
   if ( JeansMinPres )
   {
      const real Pres0 = Out[4];
      Out[4] = Hydro_CheckMinPres( Pres0, JeansMinPres_Coeff*SQR(Out[0]) );

//    recompute internal energy to be consistent with the updated pressure
      if ( EintOut != NULL  &&  Out[4] != Pres0 )
         *EintOut = EoS_DensPres2Eint( Out[0], Out[4], In+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   }


// passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass density of target passive scalars to mass fraction
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Out[ NCOMP_FLUID + FracIdx[v] ] *= _Rho;
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
//                   when FracPassive is on
//                   --> See the input parameters "FracPassive, NFrac, FracIdx"
//                3. In[] and Out[] must NOT point to the same array
//                4. In[] and Out[] should have the size of NCOMP_TOTAL_PLUS_MAG
//                5. Convert pressure to internal energy using the input EoS routine by default
//                   --> But one can also specify internal energy directly through *EintIn*,
//                       by which no EoS conversion is required and the input pressure will be useless
//                       --> Mainly used by the option LR_EINT in data reconstruction
//
// Parameter   :  In                : Array storing the input primitive variables
//                Out               : Array to store the output conserved variables
//                FracPassive       : true --> input passive scalars are mass fraction instead of density
//                NFrac             : Number of passive scalars for the option "FracPassive"
//                FracIdx           : Target variable indices for the option "FracPassive"
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensPres2Eint()
//                EoS_Table         : EoS tables for EoS_DensPres2Eint()
//                EintIn            : Pointer storing the input internal energy (see the note above)
//                                    --> Do nothing if it is NULL
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Pri2Con( const real In[], real Out[], const bool FracPassive, const int NFrac, const int FracIdx[],
                    const EoS_DP2E_t EoS_DensPres2Eint, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn )
{

   real Eint, Emag=NULL_REAL;

// passive scalars
// --> do it before invoking EoS_DensPres2Eint() since the latter requires the mass density
//     instead of mass fraction of passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// copy all passive scalars
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Out[v] = In[v];

// convert the mass fraction of target passive scalars back to mass density
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Out[ NCOMP_FLUID + FracIdx[v] ] *= In[0];
#  endif


// primitive --> conserved
   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];

#  ifdef MHD
   const real Bx = In[ MAG_OFFSET + 0 ];
   const real By = In[ MAG_OFFSET + 1 ];
   const real Bz = In[ MAG_OFFSET + 2 ];
   Emag   = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  endif
   Eint   = ( EintIn == NULL ) ? EoS_DensPres2Eint( In[0], In[4], Out+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                    EoS_AuxArray_Int, EoS_Table, NULL )
                               : *EintIn;
   Out[4] = Hydro_ConEint2Etot( Out[0], Out[1], Out[2], Out[3], Eint, Emag );


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
//                PresIn            : Pointer storing the input pressure (see the note above)
//                                    --> Do nothing if it is NULL
//
// Return      :  Flux[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real* const PresIn )
{

   const bool CheckMinPres_Yes = true;
   real InRot[ NCOMP_FLUID + NCOMP_MAG ];    // no need to include passive scalars since they don't have to be rotated

   for (int v=0; v<NCOMP_FLUID; v++)   InRot[v] = In[v];

#  ifdef MHD
   for (int v=NCOMP_FLUID; v<NCOMP_FLUID+NCOMP_MAG; v++)    InRot[v] = In[ v - NCOMP_FLUID + MAG_OFFSET ];
#  endif

   Hydro_Rotate3D( InRot, XYZ, true, NCOMP_FLUID );

#  ifdef MHD
   const real Bx   = InRot[ NCOMP_FLUID + 0 ];
   const real By   = InRot[ NCOMP_FLUID + 1 ];
   const real Bz   = InRot[ NCOMP_FLUID + 2 ];
   const real Emag = (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#  else
   const real Emag = NULL_REAL;
#  endif
   const real Pres = ( PresIn == NULL ) ? Hydro_Con2Pres( InRot[0], InRot[1], InRot[2], InRot[3], InRot[4], In+NCOMP_FLUID,
                                                          CheckMinPres_Yes, MinPres, Emag, EoS_DensEint2Pres,
                                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL )
                                        : *PresIn;
   const real _Rho = (real)1.0 / InRot[0];
   const real Vx   = _Rho*InRot[1];

   Flux[0] = InRot[1];
   Flux[1] = Vx*InRot[1] + Pres;
   Flux[2] = Vx*InRot[2];
   Flux[3] = Vx*InRot[3];
   Flux[4] = Vx*( InRot[4] + Pres );

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
// Function    :  Hydro_CheckMinPres
// Description :  Check if the input pressure is greater than the minimum allowed threshold
//
// Note        :  1. This function is used to correct unphysical (usually negative) pressure caused by
//                   numerical errors
//                   --> Usually happen in regions with high mach numbers
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
// Function    :  Hydro_CheckMinTemp
// Description :  Similar to Hydro_CheckMinPres() except that this function checks the gas temperature
//                instead of pressure
//
// Note        :  1. See Hydro_CheckMinPres()
//
// Parameter   :  InTemp  : Input temperature to be corrected
//                MinTemp : Minimum allowed temperature
//
// Return      :  InTemp != NaN --> max( InTemp, MinTemp )
//                       == NaN --> NaN
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinTemp( const real InTemp, const real MinTemp )
{

// call FMAX() only if InTemp is not NaN
   if ( InTemp == InTemp )    return FMAX( InTemp, MinTemp );
   else                       return InTemp;

} // FUNCTION : Hydro_CheckMinTemp



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
//                Emag     : Magnetic energy density (0.5*B^2) --> For MHD only
//
// Return      :  Total energy density with internal energy density greater than a given threshold
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_CheckMinEintInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real InEngy,
                               const real MinEint, const real Emag )
{

   const bool CheckMinEint_No = false;
   real InEint, OutEint, OutEngy;

   InEint  = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, InEngy, CheckMinEint_No, NULL_REAL, Emag );
   OutEint = Hydro_CheckMinEint( InEint, MinEint );

// do not modify energy (even the round-off errors) if the input data pass the check
   if ( InEint == OutEint )   OutEngy = InEngy;
   else                       OutEngy = InEngy - InEint + OutEint;

   return OutEngy;

} // FUNCTION : Hydro_CheckMinEintInEngy



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_CheckUnphysical
// Description :  Case1  :  Check if the input value is <= 0.0           (also check whether it's Inf or NAN)
//                          --> Cons and Prim must be NULL
//                Case2  :  Check if conserved variables are unphysical  (also check the existence of root in SRHD)
//                          --> Input and Prim must be NULL
//                Case3  :  Check if primitive variables are unphysical
//                          --> Input and Cons must be NULL
//
// Note        :  Two of Cons, Prim, and Input must be NULL
//                --> Can be used to check whether the values of density and pressure are unphysical
//
// Parameter   :  Cons      :  Conserved variables
//                Prim      :  Primitive variables
//                Info      :  Additional information
//                Input     :  Input value
//                Passive   :  Passive scalars
//                File      :  __FILE__
//                Function  :  __FUNCTION__
//                Line      :  __LINE__
//                Show      :  true  --> Show error message
//                             false --> Show nothing
//
// Return      :  Case1  :  true  --> Input <= 0.0  ||  >= __FLT_MAX__  ||  != itself (Nan)
//                          false --> otherwise
//                Case2  :  true  --> conserved variables are unphysical or the root(SRHD only) do not exist 
//                          false --> otherwise
//                Case3  :  true  --> primitive variables are unphysical
//                          false --> otherwise
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
bool Hydro_CheckUnphysical( const real Cons[], const real Prim[], const real* const Input, const real Passive[],
                            const char Info[], const char File[], const char Function[], const int Line, bool Show )
{

#  ifdef SRHD
   real Msqr, Dsqr, E_D, M_D, Temp, Discriminant;
#  endif

   if ( Cons == NULL && Prim == NULL && Input != NULL )
   {
      if ( *Input <= (real)0.0  ||  *Input >= __FLT_MAX__  || *Input != *Input )           goto GOTCHA;
   }

//--------------------------------------------------------------//
//------------ only check conserved variables-------------------//
//--------------------------------------------------------------//
   else if ( Cons != NULL && Prim == NULL && Input == NULL )
   {
//    check NaN
      if (  Cons[DENS] != Cons[DENS]
         || Cons[MOMX] != Cons[MOMX]
         || Cons[MOMY] != Cons[MOMY]
         || Cons[MOMZ] != Cons[MOMZ]
         || Cons[ENGY] != Cons[ENGY]  )                                                   goto GOTCHA;

//    check +inf and -inf for MOMX/Y/Z and negative for DENS and ENGY
      if (  (real)  TINY_NUMBER >= Cons[DENS] || Cons[DENS]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Cons[MOMX] || Cons[MOMX]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Cons[MOMY] || Cons[MOMY]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Cons[MOMZ] || Cons[MOMZ]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Cons[ENGY] || Cons[ENGY]  >= (real)HUGE_NUMBER )       goto GOTCHA;


//    calculate Discriminant
#     ifdef SRHD
      Msqr         = SQR(Cons[MOMX]) + SQR(Cons[MOMY]) + SQR(Cons[MOMZ]);
      Dsqr         = SQR(Cons[0]);
      E_D          = Cons[4] / Cons[0];
      M_D          = SQRT( Msqr / Dsqr );
      Temp         = SQRT( E_D*E_D + (real)2.0*E_D );
      Discriminant = ( Temp + M_D ) * ( Temp - M_D );

//    check Discriminant
      if ( Discriminant <= TINY_NUMBER )                                                  goto GOTCHA;
#     endif

//    pass all checks
      return false;
   }

//--------------------------------------------------------------//
//------------ only check primitive variables-------------------//
//--------------------------------------------------------------//

   else if ( Cons == NULL && Prim != NULL && Input == NULL )
   {
//    check NaN
      if (  Prim[0] != Prim[0]
         || Prim[1] != Prim[1]
         || Prim[2] != Prim[2]
         || Prim[3] != Prim[3]
         || Prim[4] != Prim[4]  )                                                         goto GOTCHA;

//    check +inf and -inf for velocities and negative for mass density and pressure
      if (  (real)  TINY_NUMBER >= Prim[0] || Prim[0]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Prim[1] || Prim[1]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Prim[2] || Prim[2]  >= (real)HUGE_NUMBER
         || (real) -HUGE_NUMBER >= Prim[3] || Prim[3]  >= (real)HUGE_NUMBER
         || (real)  TINY_NUMBER >= Prim[4] || Prim[4]  >= (real)HUGE_NUMBER )             goto GOTCHA;


//    pass all checks
      return false;
   }
   else
   {
      printf("ERROR: Two of Cons, Prim, and Input must be NULL !! at file <%s>, line <%d>, function <%s>\n",
              File, Line, Function );
      return true;
    }

//  print all variables if goto GOTCHA
    GOTCHA:
    {
      if ( Show )
       {

         if ( Cons == NULL && Prim == NULL && Input != NULL && Info != NULL )
         {
            printf( "ERROR: invalid %s (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                     Info, *Input, File, Line, Function );
#           if ( NCOMP_PASSIVE > 0 )
            if ( Passive != NULL )
            {
              printf( "        Passive scalars:" );
              for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] ); 
              printf( "\n" );
            }
#           endif
         }
         else if ( Cons != NULL && Prim == NULL && Input == NULL )
         {
            printf( "ERROR: unphysical conserved variables at file <%s>, line <%d>, function <%s>\n",
                            File, Line, Function );
#           ifdef SRHD
            printf( "       D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e, E^2+2*E*D-|M|^2=%14.7e\n",
                            Cons[DENS], Cons[MOMX], Cons[MOMY], Cons[MOMZ], Cons[ENGY], Discriminant );
#           else
            printf( "       D=%14.7e, Mx=%14.7e, My=%14.7e, Mz=%14.7e, E=%14.7e\n",
                            Cons[DENS], Cons[MOMX], Cons[MOMY], Cons[MOMZ], Cons[ENGY] );
#           endif
         }
         else if ( Cons == NULL && Prim != NULL && Input == NULL )
         {
            printf( "ERROR: unphysical primitive variables at file <%s>, line <%d>, function <%s>\n",
                            File, Line, Function );
            printf( "       rho=%14.7e, Ux=%14.7e, Uy=%14.7e, Uz=%14.7e, P=%14.7e\n",
                            Prim[0], Prim[1], Prim[2], Prim[3], Prim[4] );
         }

       }

      return true;
    }
}



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
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut )
{

   const bool CheckMinEint_No = false;
   real Eint, Pres;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, Emag );
   Pres = EoS_DensEint2Pres( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

   if ( CheckMinPres )   Pres = Hydro_CheckMinPres( Pres, MinPres );

   if ( EintOut != NULL )  *EintOut = Eint;

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
                     const bool CheckMinEint, const real MinEint, const real Emag )
{

//###NOTE: assuming Etot = Eint + Ekin + Emag
   real Eint;

   Eint  = Engy - (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / Dens;
#  ifdef MHD
   Eint -= Emag;
#  endif

   if ( CheckMinEint )   Eint = Hydro_CheckMinEint( Eint, MinEint );

   return Eint;

} // FUNCTION : Hydro_Con2Eint



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
// Function    :  Hydro_Con2Temp
// Description :  Evaluate the fluid temperature
//
// Note        :  1. Invoke the EoS routine EoS_DensEint2Temp() to support different EoS
//                2. Temperature is in kelvin
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Energy density
//                Passive           : Passive scalars
//                CheckMinTemp      : Apply temperature floor by calling Hydro_CheckMinTemp()
//                                    --> In some cases we actually want to check if temperature becomes unphysical,
//                                        for which we don't want to enable this option
//                MinTemp           : Temperature floor
//                Bmag              : Magnetic energy density (0.5*B^2) --> For MHD only
//                EoS_DensEint2Temp : EoS routine to compute the gas temperature
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Temp()
//                EoS_Table         : EoS tables for EoS_DensEint2Temp()
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_Con2Temp( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinTemp, const real MinTemp, const real Emag,
                     const EoS_DE2T_t EoS_DensEint2Temp, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( EoS_DensEint2Temp == NULL )
   {
#     ifdef __CUDACC__
      printf( "ERROR : EoS_DensEint2Temp == NULL at file <%s>, line <%d>, function <%s> !!\n",
              __FILE__, __LINE__, __FUNCTION__ );
#     else
      Aux_Error( ERROR_INFO, "EoS_DensEint2Temp == NULL !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


   const bool CheckMinEint_No = false;
   real Eint, Temp;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, Emag );
   Temp = EoS_DensEint2Temp( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

   if ( CheckMinTemp )   Temp = Hydro_CheckMinTemp( Temp, MinTemp );

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



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_FLUUTILITY__
