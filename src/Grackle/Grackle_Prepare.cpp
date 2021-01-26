#include "GAMER.h"

#ifdef SUPPORT_GRACKLE


// global variables for accessing h_Che_Array[]
// --> declared in Init_MemAllocate_Grackle.cpp
extern int Che_NField;
extern int CheIdx_Dens;
extern int CheIdx_sEint;
extern int CheIdx_Ent;
extern int CheIdx_e;
extern int CheIdx_HI;
extern int CheIdx_HII;
extern int CheIdx_HeI;
extern int CheIdx_HeII;
extern int CheIdx_HeIII;
extern int CheIdx_HM;
extern int CheIdx_H2I;
extern int CheIdx_H2II;
extern int CheIdx_DI;
extern int CheIdx_DII;
extern int CheIdx_HDI;
extern int CheIdx_Metal;




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Prepare
// Description :  Fill up the input host array h_Che_Array[] for the Grackle solver
//
// Note        :  1. Prepare Che_NField variables
//                   --> Che_NField and the corresponding array indices in h_Che_Array[] (e.g., CheIdx_Dens)
//                       are declared and set by Init_MemAllocate_Grackle()
//                2. This function always prepares the latest FluSg data
//
// Parameter   :  lv          : Target refinement level
//                h_Che_Array : Host array to store the prepared data
//                NPG         : Number of patch groups prepared at a time
//                PID0_List   : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Prepare( const int lv, real h_Che_Array[], const int NPG, const int *PID0_List )
{

// check
#  ifdef GAMER_DEBUG
   if ( CheIdx_Dens == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_Dens is undefined !!\n" );
   if ( CheIdx_sEint == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_sEint is undefined !!\n" );
   if ( CheIdx_Ent == Idx_Undefined )
      Aux_Error( ERROR_INFO, "CheIdx_Ent is undefined !!\n" );

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
      if (  Idx_e == Idx_Undefined  ||  CheIdx_e == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_e is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HI == Idx_Undefined  ||  CheIdx_HI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HI is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HII == Idx_Undefined  ||  CheIdx_HII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeI == Idx_Undefined  ||  CheIdx_HeI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeI is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeII == Idx_Undefined  ||  CheIdx_HeII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
      if (  Idx_HeIII == Idx_Undefined  ||  CheIdx_HeIII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HeIII is undefined for \"GRACKLE_PRI_CHE_NSPE6\" !!\n" );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
      if (  Idx_HM == Idx_Undefined  ||  CheIdx_HM == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HM is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
      if (  Idx_H2I == Idx_Undefined  ||  CheIdx_H2I == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_H2I is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
      if (  Idx_H2II == Idx_Undefined  ||  CheIdx_H2II == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_H2II is undefined for \"GRACKLE_PRI_CHE_NSPE9\" !!\n" );
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
      if (  Idx_DI == Idx_Undefined  ||  CheIdx_DI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_DI is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
      if (  Idx_DII == Idx_Undefined  ||  CheIdx_DII == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_DII is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
      if (  Idx_HDI == Idx_Undefined  ||  CheIdx_HDI == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_HDI is undefined for \"GRACKLE_PRI_CHE_NSPE12\" !!\n" );
   }

   if ( GRACKLE_METAL ) {
      if (  Idx_Metal == Idx_Undefined  ||  CheIdx_Metal == Idx_Undefined  )
         Aux_Error( ERROR_INFO, "[Che]Idx_Metal is undefined for \"GRACKLE_METAL\" !!\n" );
   }
#  endif // #ifdef GAMER_DEBUG


   const int  Size1pg          = CUBE(PS2);
   const int  Size1v           = NPG*Size1pg;
   const real MassRatio_pe    = Const_mp / Const_me;
#  ifdef DUAL_ENERGY
   const bool CheckMinPres_No  = false;
#  else
   const bool CheckMinEint_Yes = true;
#  endif

   real *Ptr_Dens0  = h_Che_Array + CheIdx_Dens *Size1v;
   real *Ptr_sEint0 = h_Che_Array + CheIdx_sEint*Size1v;
   real *Ptr_Ent0   = h_Che_Array + CheIdx_Ent  *Size1v;
   real *Ptr_e0     = h_Che_Array + CheIdx_e    *Size1v;
   real *Ptr_HI0    = h_Che_Array + CheIdx_HI   *Size1v;
   real *Ptr_HII0   = h_Che_Array + CheIdx_HII  *Size1v;
   real *Ptr_HeI0   = h_Che_Array + CheIdx_HeI  *Size1v;
   real *Ptr_HeII0  = h_Che_Array + CheIdx_HeII *Size1v;
   real *Ptr_HeIII0 = h_Che_Array + CheIdx_HeIII*Size1v;
   real *Ptr_HM0    = h_Che_Array + CheIdx_HM   *Size1v;
   real *Ptr_H2I0   = h_Che_Array + CheIdx_H2I  *Size1v;
   real *Ptr_H2II0  = h_Che_Array + CheIdx_H2II *Size1v;
   real *Ptr_DI0    = h_Che_Array + CheIdx_DI   *Size1v;
   real *Ptr_DII0   = h_Che_Array + CheIdx_DII  *Size1v;
   real *Ptr_HDI0   = h_Che_Array + CheIdx_HDI  *Size1v;
   real *Ptr_Metal0 = h_Che_Array + CheIdx_Metal*Size1v;


#  pragma omp parallel
   {

// thread-private variables
   int  idx_p, idx_pg, PID, PID0, offset;    // idx_p/idx_pg: array indices within a patch/patch group
   real Dens, Etot, Eint, Ent;
#  ifdef DUAL_ENERGY
   real Pres;
#  else
   real Px, Py, Pz, Emag=NULL_REAL;
#  endif // #ifdef DUAL_ENERGY ... else ...
   real (*fluid)[PS1][PS1][PS1]=NULL;

   real *Ptr_Dens=NULL, *Ptr_sEint=NULL, *Ptr_Ent=NULL, *Ptr_e=NULL, *Ptr_HI=NULL, *Ptr_HII=NULL;
   real *Ptr_HeI=NULL, *Ptr_HeII=NULL, *Ptr_HeIII=NULL, *Ptr_HM=NULL, *Ptr_H2I=NULL, *Ptr_H2II=NULL;
   real *Ptr_DI=NULL, *Ptr_DII=NULL, *Ptr_HDI=NULL, *Ptr_Metal=NULL;

#  pragma omp for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0      = PID0_List[TID];
      idx_pg    = 0;
      offset    = TID*Size1pg;

      Ptr_Dens  = Ptr_Dens0  + offset;
      Ptr_sEint = Ptr_sEint0 + offset;
      Ptr_Ent   = Ptr_Ent0   + offset;
      Ptr_e     = Ptr_e0     + offset;
      Ptr_HI    = Ptr_HI0    + offset;
      Ptr_HII   = Ptr_HII0   + offset;
      Ptr_HeI   = Ptr_HeI0   + offset;
      Ptr_HeII  = Ptr_HeII0  + offset;
      Ptr_HeIII = Ptr_HeIII0 + offset;
      Ptr_HM    = Ptr_HM0    + offset;
      Ptr_H2I   = Ptr_H2I0   + offset;
      Ptr_H2II  = Ptr_H2II0  + offset;
      Ptr_DI    = Ptr_DI0    + offset;
      Ptr_DII   = Ptr_DII0   + offset;
      Ptr_HDI   = Ptr_HDI0   + offset;
      Ptr_Metal = Ptr_Metal0 + offset;

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID   = PID0 + LocalID;
         idx_p = 0;
         fluid = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;

         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
         {
            Dens = *( fluid[DENS][0][0] + idx_p );
            Etot = *( fluid[ENGY][0][0] + idx_p );

//          use the dual-energy variable to calculate the internal energy if applicable
#           ifdef DUAL_ENERGY

#           if   ( DUAL_ENERGY == DE_ENPY )
            Pres = Hydro_DensEntropy2Pres( Dens, *(fluid[ENPY][0][0]+idx_p), EoS_AuxArray_Flt[1], CheckMinPres_No, NULL_REAL );
//          EOS_GAMMA does not involve passive scalars
            Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

#           else // #ifdef DUAL_ENERGY

            Px   = *( fluid[MOMX][0][0] + idx_p );
            Py   = *( fluid[MOMY][0][0] + idx_p );
            Pz   = *( fluid[MOMZ][0][0] + idx_p );
#           ifdef MHD
            Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#           endif
            Eint = Hydro_Con2Eint( Dens, Px, Py, Pz, Etot, NULL, NULL, CheckMinEint_Yes, MIN_EINT, Emag );

#           endif // #ifdef DUAL_ENERGY ... else

//          mandatory fields
            Ptr_Dens [idx_pg] = Dens;
            Ptr_sEint[idx_pg] = Eint / Dens;
            Ptr_Ent  [idx_pg] = Etot - Eint; // non-thermal energy density

//          6-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
            Ptr_e    [idx_pg] = *( fluid[Idx_e    ][0][0] + idx_p ) * MassRatio_pe;
            Ptr_HI   [idx_pg] = *( fluid[Idx_HI   ][0][0] + idx_p );
            Ptr_HII  [idx_pg] = *( fluid[Idx_HII  ][0][0] + idx_p );
            Ptr_HeI  [idx_pg] = *( fluid[Idx_HeI  ][0][0] + idx_p );
            Ptr_HeII [idx_pg] = *( fluid[Idx_HeII ][0][0] + idx_p );
            Ptr_HeIII[idx_pg] = *( fluid[Idx_HeIII][0][0] + idx_p );
            }

//          9-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
            Ptr_HM   [idx_pg] = *( fluid[Idx_HM   ][0][0] + idx_p );
            Ptr_H2I  [idx_pg] = *( fluid[Idx_H2I  ][0][0] + idx_p );
            Ptr_H2II [idx_pg] = *( fluid[Idx_H2II ][0][0] + idx_p );
            }

//          12-species network
            if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
            Ptr_DI   [idx_pg] = *( fluid[Idx_DI   ][0][0] + idx_p );
            Ptr_DII  [idx_pg] = *( fluid[Idx_DII  ][0][0] + idx_p );
            Ptr_HDI  [idx_pg] = *( fluid[Idx_HDI  ][0][0] + idx_p );
            }

//          metallicity for metal cooling
            if ( GRACKLE_METAL )
            Ptr_Metal[idx_pg] = *( fluid[Idx_Metal][0][0] + idx_p );

            idx_p  ++;
            idx_pg ++;
         } // i,j,k

      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

   } // end of OpenMP parallel region


// set cell size and link pointers for different fields
   Che_FieldData->grid_dx         = amr->dh[lv];

   Che_FieldData->density         = Ptr_Dens0;
   Che_FieldData->internal_energy = Ptr_sEint0;

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   Che_FieldData->e_density       = Ptr_e0;
   Che_FieldData->HI_density      = Ptr_HI0;
   Che_FieldData->HII_density     = Ptr_HII0;
   Che_FieldData->HeI_density     = Ptr_HeI0;
   Che_FieldData->HeII_density    = Ptr_HeII0;
   Che_FieldData->HeIII_density   = Ptr_HeIII0;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   Che_FieldData->HM_density      = Ptr_HM0;
   Che_FieldData->H2I_density     = Ptr_H2I0;
   Che_FieldData->H2II_density    = Ptr_H2II0;
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   Che_FieldData->DI_density      = Ptr_DI0;
   Che_FieldData->DII_density     = Ptr_DII0;
   Che_FieldData->HDI_density     = Ptr_HDI0;
   }

   if ( GRACKLE_METAL )
   Che_FieldData->metal_density   = Ptr_Metal0;

} // FUNCTION : Grackle_Prepare



#endif // #ifdef SUPPORT_GRACKLE
