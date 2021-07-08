#include "GAMER.h"


void Int_MinMod1D  ( const real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                           real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool OppSign0thOrder );
void Int_MinMod3D  (       real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                           real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool OppSign0thOrder );
void Int_vanLeer(          real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                           real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool OppSign0thOrder );
void Int_CQuadratic(       real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                           real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_Quadratic (       real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
		           real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_CQuartic  (       real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
	                   real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );
void Int_Quartic   (       real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
	                   real FData[], const int FSize[3], const int FStart[3], const int NComp,
                     const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );




//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolate
// Description :  Perform spatial interpolation using different schemes
//
// Note        :  Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//
// Parameter   :  CData           : Input coarse-grid array
//                CSize           : Size of CData[]
//                CStart          : (x,y,z) starting indices to perform interpolation on CData[]
//                CRange          : Number of coarse cells along each direction to perform interpolation
//                FData           : Output fine-grid array
//                FSize           : Size of FData[]
//                FStart          : (x,y,z) starting indices to store the interpolation results
//                NComp           : Number of components in the CData and FData array
//                IntScheme       : Interpolation scheme
//                                  --> currently supported schemes include
//                                      INT_MINMOD3D : MinMod-3D
//                                      INT_MINMOD1D : MinMod-1D
//                                      INT_VANLEER  : vanLeer
//                                      INT_CQUAD    : conservative quadratic
//                                      INT_QUAD     : quadratic
//                                      INT_CQUAR    : conservative quartic
//                                      INT_QUAR     : quartic
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                Monotonic       : Ensure that all interpolation results are monotonic
//                                  --> Useful for interpolating positive-definite variables, such as density, energy, ...
//                IntMonoCoeff    : min-mod coefficient for spatial interpolations
//                OppSign0thOrder : Apply 0th-order interpolation if the values to be interpolated change
//                                  signs in adjacent cells
//                                  --> See Int_MinMod1D() for details
//-------------------------------------------------------------------------------------------------------
void Interpolate( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData[], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic[],
                  const real IntMonoCoeff, const bool OppSign0thOrder )
{

// check
#  ifdef GAMER_DEBUG
   int NGhost, NSide;
   Int_Table( IntScheme, NSide, NGhost );

   for (int d=0; d<3; d++)
   {
      if ( CSize[d] < 0 )     Aux_Error( ERROR_INFO, "CSize[%d] = %d < 0 !!\n", d, CSize[d] );
      if ( FSize[d] < 0 )     Aux_Error( ERROR_INFO, "FSize[%d] = %d < 0 !!\n", d, FSize[d] );
      if ( CStart[d] < NGhost  ||  CStart[d] >= CSize[d]-NGhost )
         Aux_Error( ERROR_INFO, "incorrect CStart[%d] = %d (Min = %d, Max = %d) !!\n",
                    d, CStart[d], NGhost, CSize[d]-NGhost-1 );
      if ( FStart[d] < 0  ||  FStart[d] >= FSize[d]-1 )
         Aux_Error( ERROR_INFO, "incorrect FStart[%d] = %d (Min = %d, Max = %d) !!\n",
                    d, FStart[d], 0, FSize[d]-2 );
      if ( CStart[d]+CRange[d] >= CSize[d]-NGhost+1 )
         Aux_Error( ERROR_INFO, "incorrect CStart[%d] (%d) + CRange[%d] (%d) = %d (Max = %d) !!\n",
                    d, CStart[d], d, CRange[d], CStart[d]+CRange[d], CSize[d]-NGhost );
   }

   if ( UnwrapPhase )
   {
#     if ( MODEL == ELBDM )
      if ( IntScheme == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "unsupported phase interpolation scheme (%d) !!\n", IntScheme );
#     else
      Aux_Error( ERROR_INFO, "phase unwrapping is useful in ELBDM model only !!\n" );
#     endif
   }
#  endif // #ifdef GAMER_DEBUG


   switch ( IntScheme )
   {
      case INT_MINMOD3D :
         Int_MinMod3D  ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase,
                         OppSign0thOrder );
         break;

      case INT_MINMOD1D :
         Int_MinMod1D  ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp,
                         OppSign0thOrder );
         break;

      case INT_VANLEER :
         Int_vanLeer   ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase,
                         OppSign0thOrder );
         break;

      case INT_CQUAD :
         Int_CQuadratic( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase, Monotonic, IntMonoCoeff,
                         OppSign0thOrder );
         break;

      case INT_QUAD :
         Int_Quadratic ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase, Monotonic, IntMonoCoeff,
                         OppSign0thOrder );
         break;

      case INT_CQUAR :
         Int_CQuartic  ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase, Monotonic, IntMonoCoeff,
                         OppSign0thOrder );
         break;

      case INT_QUAR :
         Int_Quartic   ( CData, CSize, CStart, CRange, FData, FSize, FStart, NComp, UnwrapPhase, Monotonic, IntMonoCoeff,
                         OppSign0thOrder );
         break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
   } // switch ( IntScheme )

} // FUNCTION : Interpolate
