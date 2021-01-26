#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check
// Description :  Trigger the auxiliary check functions
//-------------------------------------------------------------------------------------------------------
void Aux_Check( )
{

   if ( OPT__CK_REFINE )
      for (int lv=0; lv<NLEVEL-1; lv++)   Aux_Check_Refinement( lv, "DIAGNOSIS" );

   if ( OPT__CK_PROPER_NESTING )
      for (int lv=1; lv<NLEVEL; lv++)     Aux_Check_ProperNesting( lv, "DIAGNOSIS" );

   if ( OPT__CK_CONSERVATION )            Aux_Check_Conservation( "DIAGNOSIS" );

#  if ( NCOMP_PASSIVE > 0 )
   if ( OPT__CK_NORMALIZE_PASSIVE )
      for (int lv=0; lv<NLEVEL; lv++)     Aux_Check_NormalizePassive( lv, "DIAGNOSIS" );
#  endif

   if ( OPT__CK_RESTRICT )
      for (int lv=0; lv<NLEVEL-1; lv++)   Aux_Check_Restrict( lv, "DIAGNOSIS" );

   if ( OPT__CK_FINITE )
      for (int lv=0; lv<NLEVEL; lv++)     Aux_Check_Finite( lv, "DIAGNOSIS" );

   if ( OPT__CK_PATCH_ALLOCATE )
      for (int lv=0; lv<NLEVEL; lv++)     Aux_Check_PatchAllocate( lv, "DIAGNOSIS" );

   if ( OPT__CK_FLUX_ALLOCATE )
      for (int lv=0; lv<NLEVEL-1; lv++)   Aux_Check_FluxAllocate( lv, "DIAGNOSIS" );

#  if ( MODEL == HYDRO && !defined SRHD )
   if ( OPT__CK_NEGATIVE )
      for (int lv=0; lv<NLEVEL; lv++)     Hydro_Aux_Check_Negative( lv, OPT__CK_NEGATIVE, "DIAGNOSIS" );
#  endif

#  ifdef PARTICLE
   if ( OPT__CK_PARTICLE )                Par_Aux_Check_Particle( "DIAGNOSIS" );
#  endif

#  ifdef MHD
   if ( OPT__CK_INTERFACE_B )             MHD_Aux_Check_InterfaceB( "DIAGNOSIS" );

   if ( OPT__CK_DIVERGENCE_B )            MHD_Aux_Check_DivergenceB( OPT__CK_DIVERGENCE_B==2, "DIAGNOSIS" );
#  endif

   if ( OPT__CK_MEMFREE != 0.0 )          Aux_Check_MemFree( OPT__CK_MEMFREE, "DIAGNOSIS" );

} // FUNCTION : Aux_Check
