#include "GAMER.h"

#if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Refine
// Description :  Construct patches at level "lv+1" according to the flagging result at level "lv"
//
// Note        :  1. This function will also construct buffer patches at level "lv+1" by calling Refine_Buffer()
//                2. Data of all sibling-buffer patches must be prepared in advance for creating new
//                   fine-grid patches by spatial interpolation
//                3. If LOAD_BALANCE is turned on and UseLBFunc==true, this function will invoke LB_Refine() instead
//
// Parameter   :  lv        : Target refinement level to be refined
//                UseLBFunc : Invoke the load-balance alternative functions for the grid refinement
//                            --> USELB_YES : use the load-balance alternative functions
//                                USELB_NO  : do not use the load-balance alternative functions
//-------------------------------------------------------------------------------------------------------
void Refine( const int lv, const UseLBFunc_t UseLBFunc )
{

// invoke the load-balance refine function
#  ifdef LOAD_BALANCE
   if ( UseLBFunc == USELB_YES )
   {
      LB_Refine( lv );
      return;
   }
#  endif


// check
   if ( lv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function <%s> should NOT be applied to the finest level !!\n",
                   __FUNCTION__ );
      return;
   }

// check the synchronization
   if ( NPatchTotal[lv+1] != 0 )    Mis_CompareRealValue( Time[lv], Time[lv+1], __FUNCTION__, true );


   const int  Width       = PS1*amr->scale[lv+1];  // scale of a single patch      at level "lv+1"
   const int  CFluSg      = amr->FluSg[lv  ];      // sandglass of fluid variables at level "lv"
   const int  FFluSg      = amr->FluSg[lv+1];      // sandglass of fluid variables at level "lv+1"
#  ifdef GRAVITY
   const int  CPotSg      = amr->PotSg[lv  ];      // sandglass of potential       at level "lv"
   const int  FPotSg      = amr->PotSg[lv+1];      // sandglass of potential       at level "lv+1"
   const bool UsePot      = ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT );
#  endif
#  ifdef MHD
   const int  CMagSg      = amr->MagSg[lv  ];      // sandglass of magnetic field  at level "lv"
   const int  FMagSg      = amr->MagSg[lv+1];      // sandglass of magnetic field  at level "lv+1"
#  endif

   int *Cr            = NULL;    // corner coordinates
   int *BufGrandTable = NULL;    // table recording the patch IDs of grandson buffer patches
   int *BufSonTable   = NULL;    // table recording the linking index of each buffer father patch to BufGrandTable


// parameters for spatial interpolation
   const int CRange_CC[3] = { PS1, PS1, PS1 };
   const int FSize_CC     = PS2;
   const int FStart_CC[3] = { 0, 0, 0 };
   const int FSize_CC3[3] = { FSize_CC, FSize_CC, FSize_CC };

   int NSide_Flu, CGhost_Flu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, CGhost_Flu );

   const int CSize_Flu     = PS1 + 2*CGhost_Flu;
   const int CStart_Flu[3] = { CGhost_Flu, CGhost_Flu, CGhost_Flu };
   const int CSize_Flu3[3] = { CSize_Flu, CSize_Flu, CSize_Flu };

   real Flu_CData[NCOMP_TOTAL][CSize_Flu][CSize_Flu][CSize_Flu];  // coarse-grid fluid array for interpolation
   real Flu_FData[NCOMP_TOTAL][FSize_CC ][FSize_CC ][FSize_CC ];  // fine-grid fluid array storing the interpolation result

#  ifdef GRAVITY
   int NSide_Pot, CGhost_Pot;
   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, CGhost_Pot );

   const int CSize_Pot     = PS1 + 2*CGhost_Pot;
   const int CStart_Pot[3] = { CGhost_Pot, CGhost_Pot, CGhost_Pot };

   real Pot_CData[CSize_Pot][CSize_Pot][CSize_Pot];   // coarse-grid potential array for interpolation
   real Pot_FData[FSize_CC ][FSize_CC ][FSize_CC ];   // fine-grid potential array storing the interpolation result
#  endif

#  ifdef MHD
   const int CRange_Mag[3]    = { PS1, PS1, PS1 };
   const int FSize_Mag [3][3] = {  { PS2P1, PS2,   PS2   },
                                   { PS2,   PS2P1, PS2   },
                                   { PS2,   PS2,   PS2P1 }  };
   const int FStart_Mag[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

   int NSide_Mag_Useless, CGhost_Mag;
   Int_Table( OPT__REF_MAG_INT_SCHEME, NSide_Mag_Useless, CGhost_Mag );

   const int CSize_Mag_T      = PS1 + 2*CGhost_Mag;   // coarse-grid size along the transverse (_T) / normal (_N) direction
   const int CSize_Mag_N      = PS1P1;
   const int CStart_Mag[3][3] = {  { 0,          CGhost_Mag, CGhost_Mag },
                                   { CGhost_Mag, 0,          CGhost_Mag },
                                   { CGhost_Mag, CGhost_Mag, 0          }  };
   const int CSize_Mag [3][3] = {  { CSize_Mag_N, CSize_Mag_T, CSize_Mag_T },
                                   { CSize_Mag_T, CSize_Mag_N, CSize_Mag_T },
                                   { CSize_Mag_T, CSize_Mag_T, CSize_Mag_N }  };

   real Mag_CData[NCOMP_MAG][ CSize_Mag_N*SQR(CSize_Mag_T) ];  // coarse-grid B field array for interpolation
   real Mag_FData[NCOMP_MAG][ PS2P1*SQR(PS2) ];                // fine-grid B field array storing the interpolation result

   real *Mag_FInterface_Ptr [6] = { NULL, NULL, NULL, NULL, NULL, NULL };
   real *Mag_FInterface_Data[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
   for (int s=0; s<6; s++)    Mag_FInterface_Data[s] = new real [ SQR(PS2) ];

   bool *JustRefined = new bool [ amr->num[lv] ];
   for (int PID=0; PID<amr->num[lv]; PID++)  JustRefined[PID] = false;
#  endif // #ifdef MHD


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   const int   NDer       = 0;
   const long *DerVarList = NULL;

   int BC_Face[26], BC_Face_tmp[3], FluVarIdxList[NCOMP_TOTAL];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1  &&  OPT__BC_FLU[BC_Face_tmp[2]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1  &&  OPT__BC_FLU[BC_Face_tmp[1]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1  &&  OPT__BC_FLU[BC_Face_tmp[0]] != BC_FLU_PERIODIC )   BC_Face[s] = BC_Face_tmp[0];
      else                                                                                   BC_Face[s] = NULL_INT;
   }

   for (int v=0; v<NCOMP_TOTAL; v++)   FluVarIdxList[v] = v;



// a. record the tables "BufGrandTable" and "BufFathTable"
// *** the grandson patch indices must be stored in advance before we remove all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
   if ( lv < NLEVEL-2 )
   {
      const int NBufSon = amr->NPatchComma[lv+1][27] - amr->NPatchComma[lv+1][1];
      const int NBufFa  = amr->NPatchComma[lv  ][27] - amr->NPatchComma[lv  ][1];

      BufGrandTable = new int [NBufSon];
      BufSonTable   = new int [NBufFa ];

//    initialize the table BufSonTable as -1
//#     pragma omp parallel for
      for (int t=0; t<NBufFa; t++)  BufSonTable[t] = -1;

//#     pragma omp parallel for
      for (int m=0; m<NBufSon; m+=8)
      {
//       record the grandson patch ID
         for (int n=m; n<m+8; n++)
         {
            const int BufSonPID = amr->NPatchComma[lv+1][1] + n;

            BufGrandTable[n] = amr->patch[0][lv+1][BufSonPID]->son;
         }

//       record the index of BufGrandTable array for father patches with son
         const int BufFaID = amr->patch[0][lv+1][ amr->NPatchComma[lv+1][1] + m ]->father
                             - amr->NPatchComma[lv][1];

         BufSonTable[BufFaID] = m;
      }
   } // if ( lv < NLEVEL-2 )



// b. deallocate all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
//#  pragma omp parallel for
   for (int PID=amr->NPatchComma[lv+1][1]; PID<amr->NPatchComma[lv+1][27]; PID++)
   {
      amr->patch[0][lv+1][PID]->son = -1;
      amr->pdelete( lv+1, PID, OPT__REUSE_MEMORY );
   }

//#  pragma omp parallel for
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][27]; PID++)
      amr->patch[0][lv][PID]->son = -1;



// c. check the refinement flags for all real patches at level "lv"
// ------------------------------------------------------------------------------------------------

// (c1) construct new child patches (allocate one patch group at a time)
//      --> note that we must do this BEFORE deallocating any child patch to retain high-resolution
//          B field on the boundaries of newly allocated patches
// ================================================================================================
//#  pragma omp parallel for private( ??? )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      patch_t *Pedigree = amr->patch[0][lv][PID];  // fixed to Sg=0 for the patch relation

      if ( Pedigree->flag  &&  Pedigree->son == -1 )
      {
//       (c1.1) construct relation : father -> child
         Pedigree->son = amr->num[lv+1];


//       (c1.2) allocate child patches and construct relation : child -> father
         Cr = Pedigree->corner;

         amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, true, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, true, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, true, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, true, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, true, true, true );
         amr->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, true, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, true, true, true );
         amr->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, true, true, true );

//       record the newly refined father patches
#        ifdef MHD
         JustRefined[PID] = true;
#        endif


//       (c1.3) assign data to child patches by spatial interpolation
//       (c1.3.1) fill up the central region of CData
         int i_out, j_out, k_out;

//       fluid data
         for (int v=0; v<NCOMP_TOTAL; v++)   {
         for (int k=0; k<PS1; k++)  {  k_out = k + CGhost_Flu;
         for (int j=0; j<PS1; j++)  {  j_out = j + CGhost_Flu;
         for (int i=0; i<PS1; i++)  {  i_out = i + CGhost_Flu;

            Flu_CData[v][k_out][j_out][i_out] = amr->patch[CFluSg][lv][PID]->fluid[v][k][j][i];

         }}}}

//       potential data
#        ifdef GRAVITY
         if ( UsePot )
         for (int k=0; k<PS1; k++)  {  k_out = k + CGhost_Pot;
         for (int j=0; j<PS1; j++)  {  j_out = j + CGhost_Pot;
         for (int i=0; i<PS1; i++)  {  i_out = i + CGhost_Pot;

            Pot_CData[k_out][j_out][i_out] = amr->patch[CPotSg][lv][PID]->pot[k][j][i];

         }}}
#        endif

//       magnetic field
#        ifdef MHD
         int idx_B_in, idx_B_out;

//       Bx
         idx_B_in = 0;
         for (int k=CGhost_Mag; k<CGhost_Mag+PS1; k++)  {
         for (int j=CGhost_Mag; j<CGhost_Mag+PS1; j++)  {  idx_B_out = IDX321(          0, j, k, CSize_Mag_N, CSize_Mag_T );
         for (int i=0;          i<CSize_Mag_N;    i++)  {
            Mag_CData[MAGX][ idx_B_out ++ ] = amr->patch[CMagSg][lv][PID]->magnetic[MAGX][ idx_B_in ++ ];
         }}}

//       By
         idx_B_in = 0;
         for (int k=CGhost_Mag; k<CGhost_Mag+PS1; k++)  {
         for (int j=0;          j<CSize_Mag_N;    j++)  {  idx_B_out = IDX321( CGhost_Mag, j, k, CSize_Mag_T, CSize_Mag_N );
         for (int i=CGhost_Mag; i<CGhost_Mag+PS1; i++)  {
            Mag_CData[MAGY][ idx_B_out ++ ] = amr->patch[CMagSg][lv][PID]->magnetic[MAGY][ idx_B_in ++ ];
         }}}

//       Bz
         idx_B_in = 0;
         for (int k=0;          k<CSize_Mag_N;    k++)  {
         for (int j=CGhost_Mag; j<CGhost_Mag+PS1; j++)  {  idx_B_out = IDX321( CGhost_Mag, j, k, CSize_Mag_T, CSize_Mag_T );
         for (int i=CGhost_Mag; i<CGhost_Mag+PS1; i++)  {
            Mag_CData[MAGZ][ idx_B_out ++ ] = amr->patch[CMagSg][lv][PID]->magnetic[MAGZ][ idx_B_in ++ ];
         }}}
#        endif // #ifdef MHD


//       (c1.3.2) fill up the ghost zone of CData (no interpolation is required)
         int    loop[3], offset_out[3], offset_in[3], i_in, j_in, k_in, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3];
         double xyz_flu[3];

//       calculate the corner coordinates of the coarse-grid data for the user-specified B.C.
         for (int d=0; d<3; d++)    xyz_flu[d] = Pedigree->EdgeL[d] + (0.5-CGhost_Flu)*amr->dh[lv];

//       (c1.3.2.1) prepare the fluid data
         for (int sib=0; sib<NSide_Flu; sib++)
         {
            const int SibPID = Pedigree->sibling[sib];

            for (int d=0; d<3; d++)
            {
               loop      [d] = TABLE_01( sib, 'x'+d, CGhost_Flu, PS1, CGhost_Flu );
               offset_out[d] = TABLE_01( sib, 'x'+d, 0, CGhost_Flu, CGhost_Flu+PS1 );
            }

//          (c1.3.2.1-1) if the target sibling patch exists --> just copy data from it directly
            if ( SibPID >= 0 )
            {
               for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-CGhost_Flu, 0, 0 );

               for (int v=0; v<NCOMP_TOTAL; v++)  {
               for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
               for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
               for (int i=0; i<loop[0]; i++)  {  i_out = i + offset_out[0];  i_in = i + offset_in[0];

                  Flu_CData[v][k_out][j_out][i_out] = amr->patch[CFluSg][lv][SibPID]->fluid[v][k_in][j_in][i_in];

               }}}}
            } // if ( SibPID >= 0 )


//          (c1.3.2.1-2) if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
            {
               for (int d=0; d<3; d++)
               {
                  BC_Idx_Start[d] = offset_out[d];
                  BC_Idx_End  [d] = loop[d] + BC_Idx_Start[d] - 1;
               }

               BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

#              ifdef GAMER_DEBUG
               if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
                  Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

               if ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
                  Aux_Error( ERROR_INFO, "OPT__BC_FLU == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, PID %d, sib %d, lv %d) !!\n",
                             BC_Sibling, BC_Face[BC_Sibling], SibPID, PID, sib, lv );
#              endif

               switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
               {
#                 if ( MODEL == HYDRO )
                  case BC_FLU_OUTFLOW:
                     Hydro_BoundaryCondition_Outflow   ( Flu_CData[0][0][0], BC_Face[BC_Sibling], NCOMP_TOTAL, CGhost_Flu,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End );
                  break;

                  case BC_FLU_REFLECTING:
                     Hydro_BoundaryCondition_Reflecting( Flu_CData[0][0][0], BC_Face[BC_Sibling], NCOMP_TOTAL, CGhost_Flu,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                         FluVarIdxList, NDer, DerVarList );
                  break;
#                 endif

                  case BC_FLU_USER:
                     Flu_BoundaryCondition_User        ( Flu_CData[0][0][0],                      NCOMP_TOTAL,
                                                         CSize_Flu, CSize_Flu, CSize_Flu, BC_Idx_Start, BC_Idx_End,
                                                         FluVarIdxList, Time[lv], amr->dh[lv], xyz_flu, _TOTAL, lv );
                  break;

                  default:
                     Aux_Error( ERROR_INFO, "unsupported fluid B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

               } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
            } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//          (c1.3.2.1-3) it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
            else
               Aux_Error( ERROR_INFO, "SibPID = %d (lv %d, PID %d, Sib %d) !!\n", SibPID, lv, PID, sib );

         } // for (int sib=0; sib<NSide_Flu; sib++)


//       (c1.3.2.2) prepare the potential data
#        ifdef GRAVITY
         if ( UsePot )
         for (int sib=0; sib<NSide_Pot; sib++)
         {
            const int SibPID = Pedigree->sibling[sib];

            for (int d=0; d<3; d++)
            {
               loop      [d] = TABLE_01( sib, 'x'+d, CGhost_Pot, PS1, CGhost_Pot );
               offset_out[d] = TABLE_01( sib, 'x'+d, 0, CGhost_Pot, CGhost_Pot+PS1 );
            }

//          (c1.3.2.2-1) if the target sibling patch exists --> just copy data from it directly
            if ( SibPID >= 0 )
            {
               for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-CGhost_Pot, 0, 0 );

               for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
               for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
               for (int i=0; i<loop[0]; i++)  {  i_out = i + offset_out[0];  i_in = i + offset_in[0];

                  Pot_CData[k_out][j_out][i_out] = amr->patch[CPotSg][lv][SibPID]->pot[k_in][j_in][i_in];

               }}}
            } // if ( SibPID >= 0 )


//          (c1.3.2.2-2) if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
            {
               for (int d=0; d<3; d++)
               {
                  BC_Idx_Start[d] = offset_out[d];
                  BC_Idx_End  [d] = loop[d] + BC_Idx_Start[d] - 1;
               }

               BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

#              ifdef GAMER_DEBUG
               if ( BC_Face[BC_Sibling] < 0  ||  BC_Face[BC_Sibling] > 5 )
                  Aux_Error( ERROR_INFO, "incorrect BC_Face[%d] = %d !!\n", BC_Sibling, BC_Face[BC_Sibling] );

               if ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] == BC_FLU_PERIODIC )
                  Aux_Error( ERROR_INFO, "OPT__BC_FLU == BC_FLU_PERIODIC (BC_Sibling %d, BC_Face %d, SibPID %d, PID %d, sib %d, lv %d) !!\n",
                             BC_Sibling, BC_Face[BC_Sibling], SibPID, PID, sib, lv );
#              endif

//             extrapolate potential
               Poi_BoundaryCondition_Extrapolation( Pot_CData[0][0], BC_Face[BC_Sibling], 1, CGhost_Pot,
                                                    CSize_Pot, CSize_Pot, CSize_Pot, BC_Idx_Start, BC_Idx_End );
            }


//          (c1.3.2.1-3) it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
            else
               Aux_Error( ERROR_INFO, "SibPID = %d (lv %d, PID %d, Sib %d) !!\n", SibPID, lv, PID, sib );

         } // for (int sib=0; sib<NSide_Pot; sib++)
#        endif // #ifdef GRAVITY


//       (c1.3.2.3) prepare the magnetic field
#        ifdef MHD
//       interpolation on B field only requires ghost zones along the two transverse directions
//       --> skip sib>=6 since ghost zones along the diagonal directions are not required
         for (int sib=0; sib<6; sib++)
         {
            const int SibPID = Pedigree->sibling[sib];

            for (int d=0; d<3; d++)
            {
               loop      [d] = TABLE_01( sib, 'x'+d, CGhost_Mag, PS1, CGhost_Mag );
               offset_out[d] = TABLE_01( sib, 'x'+d, 0, CGhost_Mag, CGhost_Mag+PS1 );
            }

//          (c1.3.2.3-1) if the target sibling patch exists --> just copy data from it directly
            if ( SibPID >= 0 )
            {
               for (int d=0; d<3; d++)    offset_in[d] = TABLE_01( sib, 'x'+d, PS1-CGhost_Mag, 0, 0 );

//             Bx
               if ( sib != 0  &&  sib != 1 ) // skip the normal direction
               {
                  for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
                  for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
                                                    idx_B_in  = IDX321( 0, j_in,  k_in,  PS1P1,       PS1         );
                                                    idx_B_out = IDX321( 0, j_out, k_out, CSize_Mag_N, CSize_Mag_T );
                  for (int i=0; i<PS1P1;   i++)  {

                     Mag_CData[MAGX][ idx_B_out ++ ] = amr->patch[CMagSg][lv][SibPID]->magnetic[MAGX][ idx_B_in ++ ];

                  }}}
               }

//             By
               if ( sib != 2  &&  sib != 3 ) // skip the normal direction
               {
                  for (int k=0; k<loop[2]; k++)  {  k_out = k + offset_out[2];  k_in = k + offset_in[2];
                  for (int j=0; j<PS1P1;   j++)  {  j_out = j;                  j_in = j;
                                                    idx_B_in  = IDX321( offset_in[0],  j_in,  k_in,  PS1,         PS1P1       );
                                                    idx_B_out = IDX321( offset_out[0], j_out, k_out, CSize_Mag_T, CSize_Mag_N );
                  for (int i=0; i<loop[0]; i++)  {

                     Mag_CData[MAGY][ idx_B_out ++ ] = amr->patch[CMagSg][lv][SibPID]->magnetic[MAGY][ idx_B_in ++ ];

                  }}}
               }

//             Bz
               if ( sib != 4  &&  sib != 5 ) // skip the normal direction
               {
                  for (int k=0; k<PS1P1;   k++)  {  k_out = k;                  k_in = k;
                  for (int j=0; j<loop[1]; j++)  {  j_out = j + offset_out[1];  j_in = j + offset_in[1];
                                                    idx_B_in  = IDX321( offset_in[0],  j_in,  k_in,  PS1,         PS1         );
                                                    idx_B_out = IDX321( offset_out[0], j_out, k_out, CSize_Mag_T, CSize_Mag_T );
                  for (int i=0; i<loop[0]; i++)  {

                     Mag_CData[MAGZ][ idx_B_out ++ ] = amr->patch[CMagSg][lv][SibPID]->magnetic[MAGZ][ idx_B_in ++ ];

                  }}}
               }
            } // if ( SibPID >= 0 )


//          (c1.3.2.3-2) if the target sibling patch lies outside the simulation domain --> apply the specified B.C.
            else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
            {
//             work on one component at a time since the array sizes of different components are different
               for (int v=0; v<NCOMP_MAG; v++)
               {
//                get the normal direction
                  const int norm_dir = ( v == MAGX ) ? 0 :
                                       ( v == MAGY ) ? 1 :
                                       ( v == MAGZ ) ? 2 : -1;
#                 ifdef GAMER_DEBUG
                  if ( norm_dir == -1 )   Aux_Error( ERROR_INFO, "Target face-centered variable != MAGX/Y/Z !!\n" );
#                 endif

//                only need ghost zones along the two transverse directions
                  if ( sib == norm_dir*2  ||  sib == norm_dir*2+1 )  continue;

//                set array indices --> correspond to the **cell-centered** array
                  int FC_BC_Idx_Start[3], FC_BC_Idx_End[3], FC_BC_Size[3];
                  double xyz_mag[3];   // cell-centered corner coordinates for the user-specified magnetic field B.C.
                  for (int d=0; d<3; d++)
                  {
                     if ( d == norm_dir )
                     {
                        FC_BC_Idx_Start[d] = 0;
                        FC_BC_Idx_End  [d] = CSize_Mag_N - 2;
                        FC_BC_Size     [d] = CSize_Mag_N - 1;
                        xyz_mag        [d] = Pedigree->EdgeL[d] + 0.5*amr->dh[lv];
                     }

                     else
                     {
                        FC_BC_Idx_Start[d] = offset_out[d];
                        FC_BC_Idx_End  [d] = loop[d] + FC_BC_Idx_Start[d] - 1;
                        FC_BC_Size     [d] = CSize_Mag_T;
                        xyz_mag        [d] = Pedigree->EdgeL[d] + (0.5-CGhost_Mag)*amr->dh[lv];
                     }
                  }

                  BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;
                  real *Mag_CDataPtr[NCOMP_MAG] = { Mag_CData[0], Mag_CData[1], Mag_CData[2] };

                  switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
                  {
                     case BC_FLU_OUTFLOW:
                        MHD_BoundaryCondition_Outflow   ( Mag_CDataPtr, BC_Face[BC_Sibling], 1, CGhost_Mag,
                                                          FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                          &v );
                     break;

                     case BC_FLU_REFLECTING:
                        MHD_BoundaryCondition_Reflecting( Mag_CDataPtr, BC_Face[BC_Sibling], 1, CGhost_Mag,
                                                          FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                          &v );
                     break;

                     case BC_FLU_USER:
                        MHD_BoundaryCondition_User      ( Mag_CDataPtr, BC_Face[BC_Sibling], 1,
                                                          FC_BC_Size[0], FC_BC_Size[1], FC_BC_Size[2], FC_BC_Idx_Start, FC_BC_Idx_End,
                                                          &v, Time[lv], amr->dh[lv], xyz_mag, lv );
                     break;

                     default:
                        Aux_Error( ERROR_INFO, "unsupported MHD B.C. (%d) !!\n", OPT__BC_FLU[ BC_Face[BC_Sibling] ] );

                  } // switch ( OPT__BC_FLU[ BC_Face[BC_Sibling] ] )
               } // for (int v=0; v<NCOMP_MAG; v++)
            } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )


//          (c1.3.2.3-3) it will violate the proper-nesting condition if the flagged patch is NOT surrounded by siblings
            else
               Aux_Error( ERROR_INFO, "SibPID = %d (lv %d, PID %d, Sib %d) !!\n", SibPID, lv, PID, sib );
         } // for (int sib=0; sib<6; sib++)
#        endif // #ifdef MHD


//       (c1.3.3) collect fine-grid magnetic field on the coarse-fine interfaces
#        ifdef MHD
         const int didx_in[6]     = { PS1, 0, SQR(PS1), 0, CUBE(PS1), 0 };    // x=PS1/0, y=PS1/0, z=PS1/0 faces
         const int TDir[3][2]     = { {1, 2}, {0, 2}, {0, 1} };               // transverse directions
         const int stride_in_n[3] = { PS1P1, 1, 1 };
         const int stride_in_m[3] = { PS1P1*PS1, PS1P1*PS1, PS1 };

         for (int sib=0; sib<6; sib++)
         {
//          initialize Mag_FInterface_Ptr[sib] as NULL since MHD_InterpolateBField() uses that to
//          idenitfy the coarse-fine interfaces
            Mag_FInterface_Ptr[sib] = NULL;

            const int dir    = sib/2;  // spatial direction: (0,0,1,1,2,2)
            const int SibPID = Pedigree->sibling[sib];

//          skip non-periodic boundaries
            if      ( SibPID <= SIB_OFFSET_NONPERIODIC )    continue;
            else if ( SibPID == -1 )   Aux_Error( ERROR_INFO, "SibPID == -1 (lv %d, PID %d, Sib %d) !!\n", lv, PID, sib );

//          identify the coarse-fine boundaries
//          -> skip the sibling patches that have **just** been refined
//             -> treat the corresponding interfaces as coarse-coarse instead of coarse-fine interfaces
//                so that the interpolation results do not depend on the order of patches being refined
//             -> important for bitwise reproducibility
            const int SibSonPID0 = amr->patch[0][lv][SibPID]->son;
            if ( SibSonPID0 == -1  ||  JustRefined[SibPID] )   continue;

//          link pointer to the preallocated memory for identifying coarse-fine boundaries and storing data
            Mag_FInterface_Ptr[sib] = Mag_FInterface_Data[sib];

//          loop over the 4 sibling fine patches to collect the fine-grid B field on the C-F interfaces
            for (int t=0; t<4; t++)
            {
               const int LocalID    = TABLE_03( sib, t );
               const int SibSonPID  = SibSonPID0 + LocalID;
               const int didx_out_n = TABLE_02( LocalID, 'x'+TDir[dir][0], 0, PS1 );
               const int didx_out_m = TABLE_02( LocalID, 'x'+TDir[dir][1], 0, PS1 );

               for (int m=0; m<PS1; m++)  {  idx_B_in  = m*stride_in_m[dir] + didx_in[sib];
                                             idx_B_out = ( m + didx_out_m )*PS2 + didx_out_n;
               for (int n=0; n<PS1; n++)  {

                  Mag_FInterface_Ptr[sib][idx_B_out] = amr->patch[FMagSg][lv+1][SibSonPID]->magnetic[dir][idx_B_in];

                  idx_B_in  += stride_in_n[dir];
                  idx_B_out ++;
               }}
            } // for (int t=0; t<4; t++)
         } // for (int sib=0; sib<6; sib++)
#        endif // #ifdef MHD


//       (c1.3.4) perform spatial interpolation
         const bool PhaseUnwrapping_Yes   = true;
         const bool PhaseUnwrapping_No    = false;
         const bool Monotonicity_Yes      = true;
         const bool Monotonicity_No       = false;
         const bool IntOppSign0thOrder_No = false;

//       (c1.3.4.1) determine which variables require **monotonic** interpolation
         bool Monotonicity[NCOMP_TOTAL];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
#           if ( MODEL == HYDRO )
//          we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
            /*
            if ( v == DENS  ||  v == ENGY  ||  v >= NCOMP_FLUID )
                                             Monotonicity[v] = Monotonicity_Yes;
            else                             Monotonicity[v] = Monotonicity_No;
            */
                                             Monotonicity[v] = Monotonicity_Yes;

#           elif ( MODEL == ELBDM )
            if ( v != REAL  &&  v != IMAG )  Monotonicity[v] = Monotonicity_Yes;
            else                             Monotonicity[v] = Monotonicity_No;

#           else
#           error : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION IN THIS NEW MODEL ??
#           endif // MODEL
         }

//       (c1.3.4.2) interpolation
//       (c1.3.4.2-1) fluid
#        if ( MODEL == ELBDM )
         if ( OPT__INT_PHASE )
         {
//          get the wrapped phase (store in the REAL component)
#           ifdef GAMER_DEBUG
            ELBDM_GetPhase_DebugOnly( &Flu_CData[0][0][0][0], CSize_Flu );
#           else
            for (int k=0; k<CSize_Flu; k++)
            for (int j=0; j<CSize_Flu; j++)
            for (int i=0; i<CSize_Flu; i++)
               Flu_CData[REAL][k][j][i] = ATAN2( Flu_CData[IMAG][k][j][i], Flu_CData[REAL][k][j][i] );
#           endif

//          interpolate density
            Interpolate( &Flu_CData[DENS][0][0][0], CSize_Flu3, CStart_Flu, CRange_CC, &Flu_FData[DENS][0][0][0],
                         FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, &Monotonicity_Yes,
                         IntOppSign0thOrder_No );

//          interpolate phase
            Interpolate( &Flu_CData[REAL][0][0][0], CSize_Flu3, CStart_Flu, CRange_CC, &Flu_FData[REAL][0][0][0],
                         FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_Yes, &Monotonicity_No,
                         IntOppSign0thOrder_No );
         }

         else // if ( OPT__INT_PHASE )
         {
            for (int v=0; v<NCOMP_TOTAL; v++)
            Interpolate( &Flu_CData[v][0][0][0], CSize_Flu3, CStart_Flu, CRange_CC, &Flu_FData[v][0][0][0],
                         FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, Monotonicity,
                         IntOppSign0thOrder_No );
         }

         if ( OPT__INT_PHASE )
         {
//          retrieve real and imaginary parts
            real Amp, Phase, Rho;

            for (int k=0; k<FSize_CC; k++)
            for (int j=0; j<FSize_CC; j++)
            for (int i=0; i<FSize_CC; i++)
            {
               Phase = Flu_FData[REAL][k][j][i];
               Rho   = Flu_FData[DENS][k][j][i];

//             be careful about the negative density introduced from the round-off errors
               if ( Rho < (real)0.0 )
               {
                  Flu_FData[DENS][k][j][i] = (real)0.0;
                  Rho                      = (real)0.0;
               }

               Amp                      = SQRT( Rho );
               Flu_FData[REAL][k][j][i] = Amp*COS( Phase );
               Flu_FData[IMAG][k][j][i] = Amp*SIN( Phase );
            }
         }

#        else // #if ( MODEL == ELBDM )

         for (int v=0; v<NCOMP_TOTAL; v++)
         Interpolate( &Flu_CData[v][0][0][0], CSize_Flu3, CStart_Flu, CRange_CC, &Flu_FData[v][0][0][0],
                      FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, Monotonicity,
                      INT_OPP_SIGN_0TH_ORDER );

#        endif // #if ( MODEL == ELBDM ) ... else


//       (c1.3.4.2-2) potential
#        ifdef GRAVITY
         const int CSize_Pot_Temp[3] = { CSize_Pot, CSize_Pot, CSize_Pot };

         if ( UsePot )
         Interpolate( &Pot_CData[0][0][0], CSize_Pot_Temp, CStart_Pot, CRange_CC, &Pot_FData[0][0][0],
                      FSize_CC3, FStart_CC, 1, OPT__REF_POT_INT_SCHEME, PhaseUnwrapping_No, &Monotonicity_No,
                      IntOppSign0thOrder_No );
#        endif


//       (c1.3.4.2-3) magnetic field
#        ifdef MHD
         const real *Mag_CData_Ptr[NCOMP_MAG] = { Mag_CData[MAGX], Mag_CData[MAGY], Mag_CData[MAGZ] };
               real *Mag_FData_Ptr[NCOMP_MAG] = { Mag_FData[MAGX], Mag_FData[MAGY], Mag_FData[MAGZ] };

         MHD_InterpolateBField( Mag_CData_Ptr, CSize_Mag, CStart_Mag, CRange_Mag,
                                Mag_FData_Ptr, FSize_Mag, FStart_Mag, (const real**)Mag_FInterface_Ptr,
                                OPT__REF_MAG_INT_SCHEME, Monotonicity_Yes );
#        endif


//       (c1.3.4.3) check minimum density and pressure/internal energy
//       --> note that it's unnecessary to check negative passive scalars thanks to the monotonic interpolation
//       --> but we do renormalize passive scalars here
#        if ( MODEL == HYDRO  ||  MODEL == ELBDM )
         for (int k=0; k<FSize_CC; k++)
         for (int j=0; j<FSize_CC; j++)
         for (int i=0; i<FSize_CC; i++)
         {
//          check minimum density
            const real DensOld = Flu_FData[DENS][k][j][i];

            if ( DensOld < MIN_DENS )
            {
//             rescale wave function (unnecessary if OPT__INT_PHASE if off, in which case we will rescale all wave functions later)
#              if ( MODEL == ELBDM )
               if ( OPT__INT_PHASE )
               {
                  const real Rescale = SQRT( (real)MIN_DENS / DensOld );

                  Flu_FData[REAL][k][j][i] *= Rescale;
                  Flu_FData[IMAG][k][j][i] *= Rescale;
               }
#              endif

//             apply minimum density
               Flu_FData[DENS][k][j][i] = MIN_DENS;
            }


#           if ( MODEL == HYDRO  && !defined SRHD )
//          compute magnetic energy
#           ifdef MHD
            const real Emag = MHD_GetCellCenteredBEnergy( Mag_FData[MAGX], Mag_FData[MAGY], Mag_FData[MAGZ],
                                                          PS2, PS2, PS2, i, j, k );
#           else
            const real Emag = NULL_REAL;
#           endif

//          ensure consistency between pressure, total energy density, and the dual-energy variable
//          --> here we ALWAYS use the dual-energy variable to correct the total energy density
//          --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//              the runtime parameter DUAL_ENERGY_SWITCH here
#           ifdef DUAL_ENERGY
            const bool CheckMinPres_Yes = true;
            const real UseEnpy2FixEngy  = HUGE_NUMBER;
            char dummy;    // we do not record the dual-energy status here

            Hydro_DualEnergyFix( Flu_FData[DENS][k][j][i], Flu_FData[MOMX][k][j][i], Flu_FData[MOMY][k][j][i],
                                 Flu_FData[MOMZ][k][j][i], Flu_FData[ENGY][k][j][i], Flu_FData[ENPY][k][j][i],
                                 dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2], CheckMinPres_Yes, MIN_PRES,
                                 UseEnpy2FixEngy, Emag );

#           else // #ifdef DUAL_ENERGY

//          apply internal energy floor
            Flu_FData[ENGY][k][j][i]
               = Hydro_CheckMinEintInEngy( Flu_FData[DENS][k][j][i], Flu_FData[MOMX][k][j][i], Flu_FData[MOMY][k][j][i],
                                           Flu_FData[MOMZ][k][j][i], Flu_FData[ENGY][k][j][i], MIN_EINT, Emag );
#           endif // #ifdef DUAL_ENERGY ... else ...
#           endif // #if ( MODEL == HYDRO )


//          normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            if ( OPT__NORMALIZE_PASSIVE )
            {
               real Passive[NCOMP_PASSIVE];

               for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Flu_FData[ NCOMP_FLUID + v ][k][j][i];

               Hydro_NormalizePassive( Flu_FData[DENS][k][j][i], Passive, PassiveNorm_NVar, PassiveNorm_VarIdx );

               for (int v=0; v<NCOMP_PASSIVE; v++)    Flu_FData[ NCOMP_FLUID + v ][k][j][i] = Passive[v];
            }
#           endif

         } // i,j,k
#        endif // #if ( MODEL == HYDRO  ||  MODEL == ELBDM )


//       (c1.3.5) copy data from XXX_FData[] to patch pointers
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int SonPID = amr->num[lv+1] - 8 + LocalID;

            offset_in[0] = TABLE_02( LocalID, 'x', 0, PS1 );
            offset_in[1] = TABLE_02( LocalID, 'y', 0, PS1 );
            offset_in[2] = TABLE_02( LocalID, 'z', 0, PS1 );

//          fluid data
            for (int v=0; v<NCOMP_TOTAL; v++)  {
            for (int k=0; k<PS1; k++)  {  k_in = k + offset_in[2];
            for (int j=0; j<PS1; j++)  {  j_in = j + offset_in[1];
            for (int i=0; i<PS1; i++)  {  i_in = i + offset_in[0];

               amr->patch[FFluSg][lv+1][SonPID]->fluid[v][k][j][i] = Flu_FData[v][k_in][j_in][i_in];

            }}}}

//          potential data
#           ifdef GRAVITY
            if ( UsePot )
            for (int k=0; k<PS1; k++)  {  k_in = k + offset_in[2];
            for (int j=0; j<PS1; j++)  {  j_in = j + offset_in[1];
            for (int i=0; i<PS1; i++)  {  i_in = i + offset_in[0];

               amr->patch[FPotSg][lv+1][SonPID]->pot[k][j][i] = Pot_FData[k_in][j_in][i_in];

            }}}
#           endif

//          magnetic field
#           ifdef MHD
            const int Bwidth[3][3] = { {PS1P1, PS1, PS1}, {PS1, PS1P1, PS1}, {PS1, PS1, PS1P1} };

            for (int v=0; v<NCOMP_MAG; v++)     {  idx_B_out = 0;
            for (int k=0; k<Bwidth[v][2]; k++)  {  k_in      = k + offset_in[2];
            for (int j=0; j<Bwidth[v][1]; j++)  {  j_in      = j + offset_in[1];
                                                   idx_B_in  = IDX321( offset_in[0], j_in, k_in, FSize_Mag[v][0], FSize_Mag[v][1] );
            for (int i=0; i<Bwidth[v][0]; i++)  {

               amr->patch[FMagSg][lv+1][SonPID]->magnetic[v][ idx_B_out ++ ] = Mag_FData[v][ idx_B_in ++ ];

            }}}}
#           endif

//          rescale real and imaginary parts to get the correct density in ELBDM if OPT__INT_PHASE is off
#           if ( MODEL == ELBDM )
            real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

            if ( !OPT__INT_PHASE )
            for (int k=0; k<PS1; k++)
            for (int j=0; j<PS1; j++)
            for (int i=0; i<PS1; i++)
            {
               Real      = amr->patch[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i];
               Imag      = amr->patch[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i];
               Rho_Wrong = Real*Real + Imag*Imag;
               Rho_Corr  = amr->patch[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i];

//             be careful about the negative density introduced from the round-off errors
               if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
               {
                  amr->patch[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i] = (real)0.0;
                  Rescale = (real)0.0;
               }
               else
                  Rescale = SQRT( Rho_Corr/Rho_Wrong );

               amr->patch[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i] *= Rescale;
               amr->patch[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i] *= Rescale;
            }
#           endif
         } // for (int LocalID=0; LocalID<8; LocalID++)


//       (c1.4) pass particles from father to son
#        ifdef PARTICLE
         Par_PassParticle2Son_SinglePatch( lv, PID );
#        endif
      } // if ( Pedigree->flag  &&  Pedigree->son == -1 )
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// (c2) remove unflagged child patches (deallocate one patch group at a time)
//      --> note that we must do this AFTER allocating all new patches to retain high-resolution
//          B field on the boundaries of newly allocated patches
// ================================================================================================
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      patch_t *Pedigree = amr->patch[0][lv][PID];  // fixed to Sg=0 for the patch relation

      if ( !Pedigree->flag  &&  Pedigree->son != -1 )
      {
//       (c2.0) pass particles from sons to father
#        ifdef PARTICLE
         Par_PassParticle2Father( lv, PID );
#        endif


//       (c2.1) deallocate the unflagged child patches
         const int SonPID0 = Pedigree->son;
         const int NewPID0 = SonPID0;
         const int OldPID0 = amr->num[lv+1] - 8;

         for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)     amr->pdelete( lv+1, SonPID, OPT__REUSE_MEMORY );


//       (c2.2) construct relation : father -> son
         Pedigree->son = -1;


//       (c2.3) relink the child patch pointers so that no patch indices are skipped
         if ( NewPID0 != OldPID0 )
         {
            int NewPID, OldPID, GrandPID0, FaPID;

            for (int t=0; t<8; t++)
            {
               NewPID = NewPID0 + t;
               OldPID = OldPID0 + t;

//             swap pointers between the old and new PID
//             --> works no matter OPT__REUSE_MEMORY is on or off
//             --> note that when OPT__REUSE_MEMORY is on, we don't want to set pointers of OldPID as NULL
               for (int Sg=0; Sg<2; Sg++)
                  Aux_SwapPointer( (void**)&amr->patch[Sg][lv+1][OldPID], (void**)&amr->patch[Sg][lv+1][NewPID] );

//             re-construct relation : grandson -> son
               GrandPID0 = amr->patch[0][lv+1][NewPID]->son;
               if ( GrandPID0 != -1 )
               {
                  for (int GrandPID=GrandPID0; GrandPID<GrandPID0+8; GrandPID++)
                     amr->patch[0][lv+2][GrandPID]->father = NewPID;
               }
            }

//          re-construct relation : father -> son
            FaPID = amr->patch[0][lv+1][NewPID0]->father;
            amr->patch[0][lv][FaPID]->son = NewPID0;
         } // if ( NewPID0 != OldPID0 )
      } // if ( !Pedigree->flag  &&  Pedigree->son != -1 )
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// free memory
#  ifdef MHD
   for (int s=0; s<6; s++)    delete [] Mag_FInterface_Data[s];
   delete [] JustRefined;
#  endif

// initialize the amr->NPatchComma list for the buffer patches
   for (int m=1; m<28; m++)   amr->NPatchComma[lv+1][m] = amr->num[lv+1];


// d. refine buffer patches
// ------------------------------------------------------------------------------------------------
   Refine_Buffer( lv, BufSonTable, BufGrandTable );

// deallocate tables
   if ( lv < NLEVEL-2 )
   {
      delete [] BufGrandTable;
      delete [] BufSonTable;
   }


// e. re-construct tables and sibling relations
// ------------------------------------------------------------------------------------------------

// set up the BounP_IDMap for the level just created
   Buf_RecordBoundaryPatch( lv+1 );


// construct relation : siblings
   SiblingSearch( lv+1 );


// allocate flux arrays on levels "lv" and "lv+1"
   if ( amr->WithFlux )
   {
      Flu_AllocateFluxArray( lv );

      if ( lv < TOP_LEVEL-1 )
      Flu_AllocateFluxArray( lv+1 );
   }


// allocate electric field arrays on levels "lv" and "lv+1"
#  ifdef MHD
   if ( amr->WithElectric )
   {
      MHD_AllocateElectricArray( lv );

      if ( lv < TOP_LEVEL-1 )
      MHD_AllocateElectricArray( lv+1 );
   }
#  endif


// get the IDs of patches for sending and receiving data between neighbor ranks
   Buf_RecordExchangeDataPatchID( lv+1 );


// get the total number of patches at lv+1
   Mis_GetTotalPatchNumber( lv+1 );

} // FUNCTION : Refine



#if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetPhase_DebugOnly
// Description :  Alternative function to calculate phase in the debug mode so that the functions "Refine" and
//                "LB_Refine_AllocateNewPatch" will give EXACTLY THE SAME RESULTS (even the round-off errors
//                are the same)
//
// Note        :  It is found that it is necessary to use this alternative function to calculate phase
//                in order to make runs with "SERIAL", "NO LOAD_BALANCE", and "LOAD_BALANCE" have exactly
//                the same results
//
// Parameter   :  CData : Coarse-grid array
//                CSize : Size of CData in each direction
//-------------------------------------------------------------------------------------------------------
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize )
{

   const int CSize_1v = CSize*CSize*CSize;

   real *const CData_Dens = CData + DENS*CSize_1v;
   real *const CData_Real = CData + REAL*CSize_1v;
   real *const CData_Imag = CData + IMAG*CSize_1v;

   for (int t=0; t<CSize_1v; t++)   CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

} // FUNCTION :
#endif // #if ( MODEL == ELBDM  &&  defined GAMER_DEBUG )
