#include "GAMER.h"

#ifdef LOAD_BALANCE



void PrepareCData( const int FaLv, const int FaPID, real *const FaData,
                   const int FaSg_Flu, const int FaGhost_Flu, const int NSide_Flu,
                   const int FaSg_Pot, const int FaGhost_Pot, const int NSide_Pot,
                   const int FaSg_Mag, const int FaGhost_Mag,
                   const int BC_Face[], const int FluVarIdxList[] );
void LB_Refine_AllocateBufferPatch_Sibling( const int SonLv );
static int AllocateSonPatch( const int FaLv, const int *Cr, const int PScale, const int FaPID, real *CData,
                             const int CGhost_Flu, const int NSide_Flu, const int CGhost_Pot, const int NSide_Pot, const int CGhost_Mag,
                             const int BC_Face[], const int FluVarIdxList[],
                             const real *CFB_BFieldEachRank[], const int CFB_SibRank[], long CFB_OffsetEachRank[] );
static void DeallocateSonPatch( const int FaLv, const int FaPID, const int NNew_Real0, int NewSonPID0_Real[],
                                int SwitchIdx, int &RefineS2F_Send_NPatchTotal, int *&RefineS2F_Send_PIDList );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine_AllocateNewPatch
// Description :  Allocate/deallocate patches at FaLv+1
//
// Note        :  1. This function is invoked by LB_Refine()
//                2. Home/Away : target patches at home/not at home
//                3. Input New/DelCr1D_Away[] are sorted but NewCData_Away[] is unsorted
//                4. After invoking this function, some father-buffer patches at FaLv may become useless
//                   --> Currently we do not remove these patches for the consideration of better performance
//                5. This function will also allocate new father-buffer patches at FaLv and new sibling-buffer
//                   and father-buffer patches at FaLv+1
//                   --> Reallocating father-buffer patches at FaLv+1 is necessary because all buffer patches
//                       at FaLv+1 will be deallocated before allocating/deallocating real patches at FaLv+1
//                       --> Make implementation easier
//                6. All MPI lists are NOT reconstructed here
//                7. Several alternative functions are invoked here for better performance
//                   (e.g., LB_AllocateBufferPatch_Sibling() --> LB_Refine_AllocateBufferPatch_Sibling())
//
// Parameter   :  FaLv                  : Target refinement level to be refined
//                NNew_Home             : Number of home patches at FaLv to allocate son patches
//                NewPID_Home           : Patch indices of home patches at FaLv to allocate son patches
//                NNew_Away             : Number of away patches at FaLv to allocate son patches
//                NewCr1D_Away          : Padded 1D corner of away patches at FaLv to allocate son patches
//                NewCr1D_Away_IdxTable : Index table of NewCr1D_Away[]
//                NewCData_Away         : Coarse-grid data of away patches at FaLv to allocate son patches
//                NDel_Home             : Number of home patches at FaLv to deallocate son patches
//                DelPID_Home           : Patch indices of home patches at FaLv to deallocate son patches
//                NDel_Away             : Number of away patches at FaLv to deallocate son patches
//                DelCr1D_Away          : Padded 1D corner of away patches at FaLv to deallocate son patches
//
//                PARTICLE-only parameters (call-by-reference)
//                RefineS2F_Send_NPatchTotal : Total number of patches for exchanging particles from sons to fathers
//                RefineS2F_Send_PIDList     : Patch indices for exchanging particles from sons to fathers
//
//                MHD-only parameters
//                CFB_SibRank_Home : MPI ranks of the siblings of home patches
//                CFB_SibRank_Away : MPI ranks of the siblings of away patches
//                CFB_BField       : Coarse-fine interface B field array
//                CFB_NSibEachRank : Number of siblings receiving data from each rank (including its own rank)
//-------------------------------------------------------------------------------------------------------
void LB_Refine_AllocateNewPatch( const int FaLv, int NNew_Home, int *NewPID_Home, int NNew_Away,
                                 const ulong *NewCr1D_Away, const int *NewCr1D_Away_IdxTable, real *NewCData_Away,
                                 int NDel_Home, int *DelPID_Home, int NDel_Away, ulong *DelCr1D_Away,
                                 int &RefineS2F_Send_NPatchTotal, int *&RefineS2F_Send_PIDList,
                                 const int (*CFB_SibRank_Home)[6], const int (*CFB_SibRank_Away)[6],
                                 const real *CFB_BField, const int *CFB_NSibEachRank )
{

   const int SonLv    = FaLv + 1;
   const int GraLv    = FaLv + 2;
   const int SonNReal = amr->NPatchComma[SonLv][1];
   const int SonNBuff = amr->NPatchComma[SonLv][3] - SonNReal;
   const int FaNPatch = amr->num[FaLv];


// 1. get the matching lists for the away patches
// ==========================================================================================
   int *Match_New   = new int [NNew_Away];
   int *Match_Del   = new int [NDel_Away];
   int *DelPID_Away = new int [NDel_Away];

   Mis_Matching_int( FaNPatch, amr->LB->PaddedCr1DList[FaLv], NNew_Away, NewCr1D_Away, Match_New );
   Mis_Matching_int( FaNPatch, amr->LB->PaddedCr1DList[FaLv], NDel_Away, DelCr1D_Away, Match_Del );

   for (int t=0; t<NDel_Away; t++)
   {
#     ifdef GAMER_DEBUG
      if ( Match_Del[t] == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, away patch with Cr1D %lu found no matching !!\n",
                    FaLv, DelCr1D_Away[t] );
#     endif

      DelPID_Away[t] = amr->LB->PaddedCr1DList_IdxTable[FaLv][ Match_Del[t] ];
   }



// 2. remove all buffer patches at SonLv and backup the fluid, potential, and magnetic field data
//    --> backup these buffer data so as to minimize the MPI time after grid refinement
//        --> only need to exchange the buffer data that do not exist in the old MPI lists
//        --> used by DATA_AFTER_REFINE and POT_AFTER_REFINE in LB_GetBufferData()
// ==========================================================================================
   const int MirSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int FSg_Flu    = amr->FluSg[SonLv];
   const int FSg_Flu2   = 1 - FSg_Flu;
#  ifdef GRAVITY
   const int FSg_Pot    = amr->PotSg[SonLv];
   const int FSg_Pot2   = 1 - FSg_Pot;
#  endif
#  ifdef MHD
   const int FSg_Mag    = amr->MagSg[SonLv];
   const int FSg_Mag2   = 1 - FSg_Mag;
#  endif

   int NBufBk=0, NBufBk_Dup;  // BufBk : backup the data of buffer patches
                              // must set NBufBk=0 here --> otherwise it may not be initialized if SonNBuff == 0
   ulong *PCr1D_BufBk          = new ulong [SonNBuff];
   int   *PCr1D_BufBk_IdxTable = new int   [SonNBuff];
   int   *PID_BufBk            = NULL;

// to avoid GNU warnings "non-constant array new length must be specified without parentheses around the type-id [-Wvla]"
// --> see http://stackoverflow.com/questions/4523497/typedef-fixed-length-array
   /*
   real (**flu_BufBk)[PS1][PS1][PS1] = ( OPT__REUSE_MEMORY ) ? NULL : new ( real (*[SonNBuff])[PS1][PS1][PS1] );
#  ifdef GRAVITY
   real (**pot_BufBk)[PS1][PS1]      = ( OPT__REUSE_MEMORY ) ? NULL : new ( real (*[SonNBuff])[PS1][PS1] );
#  endif
   */
   typedef real flu_type[PS1][PS1][PS1];
   real (**flu_BufBk)[PS1][PS1][PS1]    = ( OPT__REUSE_MEMORY ) ? NULL : new flu_type *[SonNBuff];
#  ifdef GRAVITY
   typedef real pot_type[PS1][PS1];
   real (**pot_BufBk)[PS1][PS1]         = ( OPT__REUSE_MEMORY ) ? NULL : new pot_type *[SonNBuff];
#  endif
#  ifdef MHD
   typedef real mag_type[ PS1P1*SQR(PS1) ];
   real (**mag_BufBk)[ PS1P1*SQR(PS1) ] = ( OPT__REUSE_MEMORY ) ? NULL : new mag_type *[SonNBuff];
#  endif

   if ( SonNBuff != 0 )
   {
//    2-1. record the PID of buffer patches that will actually receive data
//    2-1-1. allocate PID_BufBk with the maximum possible number
      NBufBk_Dup = 0;
      for (int r=0; r<MPI_NRank; r++)
      {
         NBufBk_Dup += amr->LB->RecvH_NList[SonLv][r];
#        ifdef GRAVITY
         NBufBk_Dup += amr->LB->RecvG_NList[SonLv][r];
#        endif
      }

      PID_BufBk = new int [NBufBk_Dup];

//    2-1-2. record PID
      NBufBk_Dup = 0;
      for (int r=0; r<MPI_NRank; r++)
      {
//       fluid and magnetic field
         for (int t=0; t<amr->LB->RecvH_NList[SonLv][r]; t++)  PID_BufBk[ NBufBk_Dup ++ ] = amr->LB->RecvH_IDList[SonLv][r][t];

//       potential
#        ifdef GRAVITY
         for (int t=0; t<amr->LB->RecvG_NList[SonLv][r]; t++)  PID_BufBk[ NBufBk_Dup ++ ] = amr->LB->RecvG_IDList[SonLv][r][t];
#        endif
      }

//    2-1-3. sort the PID list and remove duplicates
      Mis_Heapsort( NBufBk_Dup, PID_BufBk, NULL );

      NBufBk = ( NBufBk_Dup > 0 ) ? 1 : 0;

      for (int t=1; t<NBufBk_Dup; t++)
         if ( PID_BufBk[t] != PID_BufBk[t-1] )  PID_BufBk[ NBufBk ++ ] = PID_BufBk[t];


//    2-2. backup fluid, potential and magnetic field data
      for (int t=0; t<NBufBk; t++)
      {
//       note that this patch may have only fluid[]/magnetic[], only pot[], or both
         const int SonPID = PID_BufBk[t];

//       2-2-1. if OPT__REUSE_MEMORY is used, backup fluid[], pot[], and magnetic[] by swapping the pointers of Sg=0/1
//              so that these temporarily stored buffer patch data won't be overwritten by newly allocated real patches
//              --> also note that we need to be careful about any routine that may swap patch pointers between
//                  different PIDs (e.g., DeallocateSonPatch)
//                  --> must swap back fluid[] with FSg_Flu2, pot[] with FSg_Pot2, and magnetic[] with FSg_Mag2
//                      (otherwise these pointers will point to wrong arrays after swapping patch pointers)
//              --> ugly trick ...
         if ( OPT__REUSE_MEMORY )
         {
            Aux_SwapPointer( (void**)&amr->patch[FSg_Flu ][SonLv][SonPID]->fluid,
                             (void**)&amr->patch[FSg_Flu2][SonLv][SonPID]->fluid );

#           ifdef GRAVITY
            Aux_SwapPointer( (void**)&amr->patch[FSg_Pot ][SonLv][SonPID]->pot,
                             (void**)&amr->patch[FSg_Pot2][SonLv][SonPID]->pot );
#           endif

#           ifdef MHD
            Aux_SwapPointer( (void**)&amr->patch[FSg_Mag ][SonLv][SonPID]->magnetic,
                             (void**)&amr->patch[FSg_Mag2][SonLv][SonPID]->magnetic );
#           endif
         }

//       2-2-2. if OPT__REUSE_MEMORY is not used, backup fluid[], pot[], and magnetic[] in temporary pointers
//              --> no need to backup pot_ext[] since it's actually useless for buffer patches
         else
         {
            flu_BufBk[t] = amr->patch[FSg_Flu][SonLv][SonPID]->fluid;
            amr->patch[FSg_Flu][SonLv][SonPID]->fluid    = NULL;

#           ifdef GRAVITY
            pot_BufBk[t] = amr->patch[FSg_Pot][SonLv][SonPID]->pot;
            amr->patch[FSg_Pot][SonLv][SonPID]->pot      = NULL;
#           endif

#           ifdef MHD
            mag_BufBk[t] = amr->patch[FSg_Mag][SonLv][SonPID]->magnetic;
            amr->patch[FSg_Mag][SonLv][SonPID]->magnetic = NULL;
#           endif
         }

//       2-2-3. store the padded 1D coordinate
         PCr1D_BufBk[t] = amr->patch[0][SonLv][SonPID]->PaddedCr1D;
      } // for (int t=0; t<NBufBk; t++)

//    2-2-4. sort PCr1D_BufBk
      Mis_Heapsort( NBufBk, PCr1D_BufBk, PCr1D_BufBk_IdxTable );


//    2-3. deallocate all buffer patches
      for (int SonPID=SonNReal; SonPID<amr->NPatchComma[SonLv][3]; SonPID++)
      {
//       2-3-1. reset sibling indices to -1
         for (int s=0; s<26; s++)
         {
            const int SibPID = amr->patch[0][SonLv][SonPID]->sibling[s];

            if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = -1;
         }

//       2-3-2. free memory (or just deactivate these patches if OPT__REUSE_MEMORY is on)
         amr->patch[0][SonLv][SonPID]->son = -1;
         amr->pdelete( SonLv, SonPID, OPT__REUSE_MEMORY );
      }

//    2-3-3. reset NPatchComma
      for (int m=2; m<28; m++)   amr->NPatchComma[SonLv][m] = SonNReal;

//    2-3-4. recalculate PaddedCr1DList
      amr->LB->PaddedCr1DList         [SonLv] = (ulong*)realloc( amr->LB->PaddedCr1DList         [SonLv],
                                                                 SonNReal*sizeof(ulong) );
      amr->LB->PaddedCr1DList_IdxTable[SonLv] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[SonLv],
                                                                 SonNReal*sizeof(int ) );

      for (int SonPID=0; SonPID<SonNReal; SonPID++)
         amr->LB->PaddedCr1DList[SonLv][SonPID] = amr->patch[0][SonLv][SonPID]->PaddedCr1D;

      Mis_Heapsort( SonNReal, amr->LB->PaddedCr1DList[SonLv], amr->LB->PaddedCr1DList_IdxTable[SonLv] );

   } // if ( SonNBuff != 0 )



// 3. allocate new real patches at SonLv
// ==========================================================================================
// 3.0 preparation
   const int Padded              = 1<<NLEVEL;
   const int BoxNScale_Padded[3] = { amr->BoxScale[0]/PS1 + 2*Padded,
                                     amr->BoxScale[1]/PS1 + 2*Padded,
                                     amr->BoxScale[2]/PS1 + 2*Padded }; //normalized and padded BoxScale
   const int PScale              = PS1*amr->scale[SonLv];  // scale of a single patch at SonLv
   const int NNew_Real0          = NNew_Home + NNew_Away;

   int FaPID, SonPID0, Cr3D[3], *Cr3D_Ptr=NULL, NNoFa=0;
   int *NewSonPID0_NoFa = new int [ NNew_Away ];         // NNew_Away is the maximum number this array can have
   int *NewSonPID0_All  = (int*)malloc( NNew_Real0*sizeof(int) );
   int *NewSonPID0_Real = NewSonPID0_All;
   int *NewSonPID0_Away = NewSonPID0_All + NNew_Home;


// parameters for spatial interpolation
   int CSize_Tot = 0;
   int NSide_Flu, CGhost_Flu, CSize_Flu;

   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, CGhost_Flu );

   CSize_Flu  = PS1 + 2*CGhost_Flu;
   CSize_Tot += NCOMP_TOTAL*CUBE( CSize_Flu );

#  ifdef GRAVITY
   int NSide_Pot, CGhost_Pot, CSize_Pot;

   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, CGhost_Pot );

   CSize_Pot  = PS1 + 2*CGhost_Pot;
   CSize_Tot += CUBE( CSize_Pot );
#  else
   const int NSide_Pot=NULL_INT, CGhost_Pot=NULL_INT;
#  endif

#  ifdef MHD
   int NSide_Mag_Useless, CGhost_Mag, CSize_Mag_T, CSize_Mag_N;

   Int_Table( OPT__REF_MAG_INT_SCHEME, NSide_Mag_Useless, CGhost_Mag );

   CSize_Mag_T = PS1 + 2*CGhost_Mag;   // coarse-grid size along the transverse (_T) / normal (_N) direction
   CSize_Mag_N = PS1P1;
   CSize_Tot   += NCOMP_MAG*CSize_Mag_N*SQR( CSize_Mag_T );
#  else
   const int CGhost_Mag=NULL_INT;
#  endif


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
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


// set parameters related to the coarse-fine interface B field
#  ifdef MHD
   const real *CFB_BFieldEachRank[MPI_NRank];
         long  CFB_OffsetEachRank[MPI_NRank];

   CFB_BFieldEachRank[0] = CFB_BField;
   for (int r=1; r<MPI_NRank; r++)
   CFB_BFieldEachRank[r] = CFB_BFieldEachRank[r-1] + CFB_NSibEachRank[r-1]*SQR( PS2 );

   for (int r=0; r<MPI_NRank; r++)  CFB_OffsetEachRank[r] = 0;

#  else
   const real **CFB_BFieldEachRank = NULL;
         long  *CFB_OffsetEachRank = NULL;
#  endif


// 3.1 home patches
   for (int t=0; t<NNew_Home; t++)
   {
      FaPID    = NewPID_Home[t];
      Cr3D_Ptr = amr->patch[0][FaLv][FaPID]->corner;

      NewSonPID0_All[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID,
                                            NULL,
                                            CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot, CGhost_Mag,
                                            BC_Face, FluVarIdxList,
                                            CFB_BFieldEachRank, CFB_SibRank_Home[t], CFB_OffsetEachRank );
   }


// 3.2 away patches
   for (int t=0; t<NNew_Away; t++)
   {
//    3.2.1 away patches without father patch
      if ( Match_New[t] == -1 )
      {
         FaPID = -1;
         Mis_Idx1D2Idx3D( BoxNScale_Padded, NewCr1D_Away[t], Cr3D );
         for (int d=0; d<3; d++)    Cr3D[d] = ( Cr3D[d] - Padded )*PS1;

         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D, PScale, FaPID,
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot, CGhost_Mag,
                                                NULL, NULL,
                                                CFB_BFieldEachRank, CFB_SibRank_Away[t], CFB_OffsetEachRank );

//       record the SonPID (with LocalID == 0 ) with no father at home
#        ifdef GAMER_DEBUG
         if ( NNoFa >= NNew_Away )
            Aux_Error( ERROR_INFO, "FaLv %d, NNoFa (%d) exceeds the maximum number (%d) !!\n",
                       FaLv, NNoFa, NNew_Away );
#        endif

         NewSonPID0_NoFa[ NNoFa ++ ] = NewSonPID0_Away[t];
      } // if ( Match_New[t] == -1 )


//    3.2.2 away patches with father patch
      else
      {
         FaPID    = amr->LB->PaddedCr1DList_IdxTable[FaLv][ Match_New[t] ];
         Cr3D_Ptr = amr->patch[0][FaLv][FaPID]->corner;

         NewSonPID0_Away[t] = AllocateSonPatch( FaLv, Cr3D_Ptr, PScale, FaPID,
                                                NewCData_Away+NewCr1D_Away_IdxTable[t]*CSize_Tot,
                                                CGhost_Flu, NSide_Flu, CGhost_Pot, NSide_Pot, CGhost_Mag,
                                                NULL, NULL,
                                                CFB_BFieldEachRank, CFB_SibRank_Away[t], CFB_OffsetEachRank );
      } // if ( Match_New[t] == -1 ) ... else ...
   } // for (int t=0; t<NNew_Away; t++)



// 4. allocate new father-buffer patches at FaLv and construct the relation son->father
// ==========================================================================================
   int  NNewFaBuf0;
   int *NewFaBufPID0 = NULL;

   LB_AllocateBufferPatch_Father( SonLv, false, NNew_Real0, NewSonPID0_Real, true, &NNewFaBuf0, &NewFaBufPID0 );

   LB_FindFather( SonLv, false, NNoFa, NewSonPID0_NoFa, false );



// 5. deallocate unflagged real patches at SonLv
// ==========================================================================================
   int SwitchIdx = NNew_Real0 - 1;  // Element in NewSonPID0_Real to be reset

// 5.1 home patches
   for (int t=0; t<NDel_Home; t++)
   {
      FaPID = DelPID_Home[t];

      DeallocateSonPatch( FaLv, FaPID, NNew_Real0, NewSonPID0_Real, SwitchIdx--,
                          RefineS2F_Send_NPatchTotal, RefineS2F_Send_PIDList );
   }

// 5.2 away patches
   for (int t=0; t<NDel_Away; t++)
   {
      FaPID = DelPID_Away[t];

      DeallocateSonPatch( FaLv, FaPID, NNew_Real0, NewSonPID0_Real, SwitchIdx--,
                          RefineS2F_Send_NPatchTotal, RefineS2F_Send_PIDList );
   }



// 6. record NPatchComma, NPatchTotal and LB_IdxList_Real at SonLv
// ==========================================================================================
   const int SonNReal_New = amr->NPatchComma[SonLv][1];

// 6.1 record the number of real patches at SonLv for each rank
   if ( SonNReal_New != amr->num[SonLv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d)\" !!\n",
                 SonLv, amr->NPatchComma[SonLv][1], SonLv, amr->num[SonLv] );

   for (int m=2; m<28; m++)   amr->NPatchComma[SonLv][m] = SonNReal_New;


// 6.2 record the total number of real patches at SonLv among all ranks
   Mis_GetTotalPatchNumber( SonLv );


// 6.3 record LB_IdxList_Real at SonLv
   if ( amr->LB->IdxList_Real         [SonLv] != NULL )  delete [] amr->LB->IdxList_Real         [SonLv];
   if ( amr->LB->IdxList_Real_IdxTable[SonLv] != NULL )  delete [] amr->LB->IdxList_Real_IdxTable[SonLv];

   amr->LB->IdxList_Real         [SonLv] = new long [SonNReal_New];
   amr->LB->IdxList_Real_IdxTable[SonLv] = new int  [SonNReal_New];

   for (int SonPID=0; SonPID<SonNReal_New; SonPID++)
      amr->LB->IdxList_Real[SonLv][SonPID] = amr->patch[0][SonLv][SonPID]->LB_Idx;

   Mis_Heapsort( SonNReal_New, amr->LB->IdxList_Real[SonLv], amr->LB->IdxList_Real_IdxTable[SonLv] );


// 6.4 check : no duplicate patches at FaLv and SonLv
#  ifdef GAMER_DEBUG
   for (int lv=FaLv; lv<=SonLv; lv++)
   {
      ulong *TempPaddedCr1D          = new ulong [ amr->num[lv] ];
      int   *TempPaddedCr1D_IdxTable = new  int  [ amr->num[lv] ];

      for (int PID=0; PID<amr->num[lv]; PID++)
         TempPaddedCr1D[PID] = amr->patch[0][lv][PID]->PaddedCr1D;

      Mis_Heapsort( amr->num[lv], TempPaddedCr1D, TempPaddedCr1D_IdxTable );

      for (int t=1; t<amr->num[lv]; t++)
      {
         if ( TempPaddedCr1D[t] == TempPaddedCr1D[t-1] )
            Aux_Error( ERROR_INFO, "duplicate patches at lv %d, PaddedCr1D %lu, PID = %d and %d !!\n",
                       lv, TempPaddedCr1D[t], TempPaddedCr1D_IdxTable[t], TempPaddedCr1D_IdxTable[t-1] );
      }

      delete [] TempPaddedCr1D;
      delete [] TempPaddedCr1D_IdxTable;
   }
#  endif



// 7. construct new father->son relation between SonLv and FaLv and new sibling relation at FaLv
// ==========================================================================================
   LB_FindSonNotHome( FaLv, true, NULL_INT, NULL );

   LB_SiblingSearch( FaLv, false, NNewFaBuf0, NewFaBufPID0 );



// 8. allocate buffer patches at SonLv
// ==========================================================================================
// 8.1 sibling-buffer patches (invoke the alternative version which should be faster)
   LB_Refine_AllocateBufferPatch_Sibling( SonLv );
// LB_AllocateBufferPatch_Sibling( SonLv );

// 8.2 father-buffer patches
   if ( GraLv < NLEVEL )
   LB_AllocateBufferPatch_Father( GraLv, true, NULL_INT, NULL, false, NULL, NULL );



// 9. construct new father<->son relation between GraLv and SonLv and new sibling relation at SonLv
// ==========================================================================================
   const int NNew_Buf0 = ( amr->NPatchComma[SonLv][3] - amr->NPatchComma[SonLv][1] ) / 8;
   const int NNew_All0 = NNew_Real0 + NNew_Buf0;
   int *NewSonPID0_Buf, *NewSonPID_All, Count;

// 9.1 record the indices of all newly-allocated buffer patches
   NewSonPID_All  = new int [ 8*NNew_All0 ];
   NewSonPID0_All = (int*)realloc( NewSonPID0_All, NNew_All0*sizeof(int) );
   NewSonPID0_Buf = NewSonPID0_All + NNew_Real0;

   for (SonPID0=amr->NPatchComma[SonLv][1], Count=0; SonPID0<amr->NPatchComma[SonLv][3]; SonPID0+=8, Count++)
      NewSonPID0_Buf[Count] = SonPID0;

// following code seems to trigger a bug in the Intel compiler version 2016.2.181
// --> replace with the simplified version below
   /*
// original version
   for (int t=0, m=0; t<NNew_All0; t++, m+=8)
   {
      SonPID0 = NewSonPID0_All[t];

      for (int SonPID=SonPID0, n=m; SonPID<SonPID0+8; SonPID++, n++)    NewSonPID_All[n] = SonPID;
   }
   */

// simplified version
   for (int t=0; t<NNew_All0; t++)
   {
      SonPID0 = NewSonPID0_All[t];

      for (int LocalID=0; LocalID<8; LocalID++)    NewSonPID_All[ t*8 + LocalID ] = SonPID0 + LocalID;
   }

// 9.2 sibling relation
   LB_SiblingSearch( SonLv, false, NNew_All0, NewSonPID0_All );

// 9.3 father <-> son relation
   if ( GraLv < NLEVEL )
   {
      LB_FindFather( GraLv, true, NULL_INT, NULL, false );

      LB_FindSonNotHome( SonLv, false, 8*NNew_All0, NewSonPID_All );
   }



// 10. restore fluid[], pot[], and magnetic[] in the buffer patches
// ==========================================================================================
   int *Match_BufBk = new int [NBufBk];
   int  MPID;

// 10.1 get the match lists
   Mis_Matching_int( amr->num[SonLv], amr->LB->PaddedCr1DList[SonLv], NBufBk, PCr1D_BufBk, Match_BufBk );

// 10.2 reset array pointers
   for (int t=0; t<NBufBk; t++)
   {
      if ( Match_BufBk[t] != -1 )
      {
         MPID = amr->LB->PaddedCr1DList_IdxTable[SonLv][ Match_BufBk[t] ];

#        ifdef GAMER_DEBUG
         if ( MPID < amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "Match_PID = %d matches to a real patch (Match[%d] = %d, SonNReal = %d) !!\n",
                       MPID, t, Match_BufBk[t], amr->NPatchComma[SonLv][1] );
#        endif

         if ( OPT__REUSE_MEMORY )
         {
            const int OldBufPID = PID_BufBk[ PCr1D_BufBk_IdxTable[t] ];

//          note that (1) we must swap poniters even if MPID == OldBufPID (because they have different Sg)
//                    (2) we store the previous buffer data in FSg_Flu2, FSg_Pot2 and FSg_Mag2 instead of FSg_Flu, FSg_Pot, and FSg_Mag
//                    (3) it's OK to leave FSg_Flu2, FSg_Pot2, and FSg_Mag2 of MPID unmodified after swapping pointers (which can thus be NULL)
//                        since they will be allocated in LB_RecordExchangeDataPatchID() if necessary
            Aux_SwapPointer( (void**)&amr->patch[FSg_Flu ][SonLv][     MPID]->fluid,
                             (void**)&amr->patch[FSg_Flu2][SonLv][OldBufPID]->fluid );

#           ifdef GRAVITY
            Aux_SwapPointer( (void**)&amr->patch[FSg_Pot ][SonLv][     MPID]->pot,
                             (void**)&amr->patch[FSg_Pot2][SonLv][OldBufPID]->pot );
#           endif

#           ifdef MHD
            Aux_SwapPointer( (void**)&amr->patch[FSg_Mag ][SonLv][     MPID]->magnetic,
                             (void**)&amr->patch[FSg_Mag2][SonLv][OldBufPID]->magnetic );
#           endif

//          we must reallocate memory for OldBufPID if it becomes a new real patch (which always have data array allocated)
//          --> for buffer patches their data arrays are not allocated by default so we don't have to do anything here
            if ( OldBufPID < amr->NPatchComma[SonLv][1] )
            {
               amr->patch[FSg_Flu2][SonLv][OldBufPID]->hnew();
#              ifdef GRAVITY
               amr->patch[FSg_Pot2][SonLv][OldBufPID]->gnew();
#              endif
#              ifdef MHD
               amr->patch[FSg_Mag2][SonLv][OldBufPID]->mnew();
#              endif
            }
         } // if ( OPT__REUSE_MEMORY )

         else
         {
//          note that it's OK to leave FSg_Flu2, FSg_Pot2, FSg_Mag2 unmodified (which can thus be NULL) since
//          it will be allocated in LB_RecordExchangeDataPatchID if necessary
            real (*flu_ptr)[PS1][PS1][PS1] = flu_BufBk[ PCr1D_BufBk_IdxTable[t] ];
            if ( flu_ptr != NULL )
               amr->patch[FSg_Flu][SonLv][MPID]->fluid = flu_ptr;

#           ifdef GRAVITY
//          don't worry about pot_ext since it's actually useless for buffer patches
//          --> after the following operation, some buffer patches may have pot != NULL but pot_ext == NULL (for FSg_Pot)
            real (*pot_ptr)[PS1][PS1] = pot_BufBk[ PCr1D_BufBk_IdxTable[t] ];
            if ( pot_ptr != NULL )
               amr->patch[FSg_Pot][SonLv][MPID]->pot = pot_ptr;
#           endif

#           ifdef MHD
            real (*mag_ptr)[ PS1P1*SQR(PS1) ] = mag_BufBk[ PCr1D_BufBk_IdxTable[t] ];
            if ( mag_ptr != NULL )
               amr->patch[FSg_Mag][SonLv][MPID]->magnetic = mag_ptr;
#           endif
         } // if ( OPT__REUSE_MEMORY ) ... else ...
      } // if ( Match_BufBk[t] != -1 )

      else if ( ! OPT__REUSE_MEMORY )
      {
         delete [] flu_BufBk[ PCr1D_BufBk_IdxTable[t] ];
#        ifdef GRAVITY
         delete [] pot_BufBk[ PCr1D_BufBk_IdxTable[t] ];
#        endif
#        ifdef MHD
         delete [] mag_BufBk[ PCr1D_BufBk_IdxTable[t] ];
#        endif
      } // if ( Match_BufBk[t] != -1 ) ... else if ...
   } // for (int t=0; t<NBufBk; t++)



// free memory
   free( NewSonPID0_All );
   delete [] Match_New;
   delete [] Match_Del;
   delete [] Match_BufBk;
   delete [] DelPID_Away;
   delete [] NewSonPID0_NoFa;
   delete [] NewSonPID_All;
   if ( NewFaBufPID0 != NULL )   delete [] NewFaBufPID0;
   delete [] PCr1D_BufBk;
   delete [] PCr1D_BufBk_IdxTable;
   delete [] PID_BufBk;
   delete [] flu_BufBk;
#  ifdef GRAVITY
   delete [] pot_BufBk;
#  endif
#  ifdef MHD
   delete [] mag_BufBk;
#  endif

} // FUNCTION : LB_Refine_AllocateNewPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  AllocateSonPatch
// Description :  Allocate eight son patches at FaLv+1
//
// Note        :  Just to avoid duplicate code segment
//
// Parameter   :  FaLv          : Target refinement level to be refined
//                Cr            : Corner coordinates of the son patch with LocalID == 0
//                PScale        : Scale of one patch at SonLv
//                FaPID         : Father patch index (can be -1 for the away patches)
//                CData         : Coarse-grid data for assigning data to son patches by spatial interpolation
//                                (initialize as NULL if father patch is home --> prepare CData here)
//                CGhost_Flu    : Ghost size of the fluid data
//                NSide_Flu     : Number of sibling directions to prepare the ghost-zone data (6/26) for the fluid data
//                CGhost_Pot    : Ghost size of the potential data
//                NSide_Pot     : Number of sibling directions to prepare the ghost-zone data (6/26) for the potential data
//                CGhost_Mag    : Ghost size of the magnetic field data
//                BC_Face       : Corresponding boundary faces (0~5) along 26 sibling directions -> for non-periodic B.C. only
//                FluVarIdxList : List of target fluid variable indices                          -> for non-periodic B.C. only
//
//                MHD-only parameters
//                CFB_BFieldEachRank : Coarse-fine interface B field array
//                CFB_SibRank        : MPI ranks of the target sibling patches
//                CFB_OffsetEachRank : Array offset of CFB_BFieldEachRank[] for each rank
//
// Return      :  SonPID with LocalID == 0, CFB_OffsetEachRank
//-------------------------------------------------------------------------------------------------------
int AllocateSonPatch( const int FaLv, const int *Cr, const int PScale, const int FaPID, real *CData,
                      const int CGhost_Flu, const int NSide_Flu, const int CGhost_Pot, const int NSide_Pot, const int CGhost_Mag,
                      const int BC_Face[], const int FluVarIdxList[],
                      const real *CFB_BFieldEachRank[], const int CFB_SibRank[], long CFB_OffsetEachRank[] )
{

   const int SonLv   = FaLv + 1;
   const int SonPID0 = amr->num[SonLv];
   bool FaIsHome = false;

// 0. check : target father patch has no son
#  ifdef GAMER_DEBUG
   if ( FaPID != -1  &&  amr->patch[0][FaLv][FaPID]->son != -1 )
      Aux_Error( ERROR_INFO, "FaLv %d, FaPID (%d) already has sons (SonPID = %d, duplicate SonPID = %d) !!\n",
                 FaLv, FaPID, amr->patch[0][FaLv][FaPID]->son, SonPID0 );
#  endif


// 1. construct relation : father -> child
   if ( FaPID != -1 )   amr->patch[0][FaLv][FaPID]->son = SonPID0;


// 2. allocate child patches and construct relation : child -> father
   amr->pnew( SonLv, Cr[0],        Cr[1],        Cr[2],        FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1],        Cr[2],        FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1]+PScale, Cr[2],        FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1],        Cr[2]+PScale, FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1]+PScale, Cr[2],        FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0],        Cr[1]+PScale, Cr[2]+PScale, FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1],        Cr[2]+PScale, FaPID, true, true, true );
   amr->pnew( SonLv, Cr[0]+PScale, Cr[1]+PScale, Cr[2]+PScale, FaPID, true, true, true );

   amr->NPatchComma[SonLv][1] += 8;


// 3. assign data to child patches by spatial interpolation
// 3.1 prepare the coarse-grid data
   int CSize_Tot = 0;

// fluid
   const int CSg_Flu   = amr->FluSg[ FaLv];
   const int FSg_Flu   = amr->FluSg[SonLv];
   const int CSize_Flu = PS1 + 2*CGhost_Flu;
   CSize_Tot += NCOMP_TOTAL*CUBE( CSize_Flu );

// potential
#  ifdef GRAVITY
   const int CSg_Pot   = amr->PotSg[ FaLv];
   const int FSg_Pot   = amr->PotSg[SonLv];
   const int CSize_Pot = PS1 + 2*CGhost_Pot;
   CSize_Tot += CUBE( CSize_Pot );
#  else
   const int CSg_Pot   = NULL_INT;
#  endif

// magnetic field
#  ifdef MHD
   const int CSg_Mag     = amr->MagSg[ FaLv];
   const int FSg_Mag     = amr->MagSg[SonLv];
   const int CSize_Mag_T = PS1 + 2*CGhost_Mag;  // coarse-grid size along the transverse (_T) / normal (_N) direction
   const int CSize_Mag_N = PS1P1;
   CSize_Tot += NCOMP_MAG*CSize_Mag_N*SQR( CSize_Mag_T );
#  else
   const int CSg_Mag     = NULL_INT;
#  endif

// "CData == NULL" if the real father patch is home
   if ( CData == NULL )
   {
#     ifdef GAMER_DEBUG
      if ( FaPID >= amr->NPatchComma[FaLv][1]  ||  FaPID == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID (%d) is NOT home (FaNReal = %d) !!\n",
                    FaLv, FaPID, amr->NPatchComma[FaLv][1] );
#     endif

      FaIsHome = true;
      CData    = new real [CSize_Tot];

      PrepareCData( FaLv, FaPID, CData,
                    CSg_Flu, CGhost_Flu, NSide_Flu, CSg_Pot, CGhost_Pot, NSide_Pot, CSg_Mag, CGhost_Mag,
                    BC_Face, FluVarIdxList );
   }


// 3.2 perform spatial interpolation
// 3.2.1 determine which variables require **monotonic** interpolation
   const bool PhaseUnwrapping_Yes   = true;
   const bool PhaseUnwrapping_No    = false;
   const bool Monotonicity_Yes      = true;
   const bool Monotonicity_No       = false;
   const bool IntOppSign0thOrder_No = false;

   bool Monotonicity[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
#     if ( MODEL == HYDRO )
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
      /*
      if ( v == DENS  ||  v == ENGY  ||  v >= NCOMP_FLUID )
                                       Monotonicity[v] = Monotonicity_Yes;
      else                             Monotonicity[v] = Monotonicity_No;
      */
                                       Monotonicity[v] = Monotonicity_Yes;

#     elif ( MODEL == ELBDM )
      if ( v != REAL  &&  v != IMAG )  Monotonicity[v] = Monotonicity_Yes;
      else                             Monotonicity[v] = Monotonicity_No;

#     else
#     error : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION IN THIS NEW MODEL ??
#     endif // MODEL
   }

// 3.2.2 interpolation
// 3.2.2-1. fluid
   const int FSize_CC      = PS2;
   const int FSize_CC3 [3] = { FSize_CC, FSize_CC, FSize_CC };
   const int FStart_CC [3] = { 0, 0, 0 };
   const int CSize_Flu3[3] = { CSize_Flu, CSize_Flu, CSize_Flu };
   const int CStart_Flu[3] = { CGhost_Flu, CGhost_Flu, CGhost_Flu };
   const int CRange_CC [3] = { PS1, PS1, PS1 };
   const int CSize_Flu1v   = CUBE( CSize_Flu );

   real *CData_Next       = CData;
   real *const CData_Flu  = CData_Next;
#  if ( MODEL == ELBDM )
   real *const CData_Dens = CData_Flu + DENS*CSize_Flu1v;
   real *const CData_Real = CData_Flu + REAL*CSize_Flu1v;
   real *const CData_Imag = CData_Flu + IMAG*CSize_Flu1v;
#  endif
   CData_Next += NCOMP_TOTAL*CSize_Flu1v;

   real (*FData_Flu)[FSize_CC][FSize_CC][FSize_CC] = new real [NCOMP_TOTAL][FSize_CC][FSize_CC][FSize_CC];

#  if ( MODEL == ELBDM )
   if ( OPT__INT_PHASE )
   {
//    get the wrapped phase (store in the REAL component)
      for (int t=0; t<CSize_Flu1v; t++)   CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

//    interpolate density
      Interpolate( CData_Dens, CSize_Flu3, CStart_Flu, CRange_CC, &FData_Flu[DENS][0][0][0],
                   FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, &Monotonicity_Yes,
                   IntOppSign0thOrder_No );

//    interpolate phase
      Interpolate( CData_Real, CSize_Flu3, CStart_Flu, CRange_CC, &FData_Flu[REAL][0][0][0],
                   FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_Yes, &Monotonicity_No,
                   IntOppSign0thOrder_No );
   }

   else // if ( OPT__INT_PHASE )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      Interpolate( CData_Flu+v*CSize_Flu1v, CSize_Flu3, CStart_Flu, CRange_CC, &FData_Flu[v][0][0][0],
                   FSize_CC3, FStart_CC, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, Monotonicity,
                   IntOppSign0thOrder_No );
   }

   if ( OPT__INT_PHASE )
   {
//    retrieve real and imaginary parts
      real Amp, Phase, Rho;

      for (int k=0; k<FSize_CC; k++)
      for (int j=0; j<FSize_CC; j++)
      for (int i=0; i<FSize_CC; i++)
      {
         Phase = FData_Flu[REAL][k][j][i];
         Rho   = FData_Flu[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho < (real)0.0 )
         {
            FData_Flu[DENS][k][j][i] = (real)0.0;
            Rho                      = (real)0.0;
         }

         Amp                      = SQRT( Rho );
         FData_Flu[REAL][k][j][i] = Amp*COS( Phase );
         FData_Flu[IMAG][k][j][i] = Amp*SIN( Phase );
      }
   }

#  else // #if ( MODEL == ELBDM )

   InterpolateWithAdaptiveMinmod( CData_Flu, CSize_Flu3, CStart_Flu, CRange_CC, &FData_Flu[0][0][0][0],
                                  FSize_CC3, FStart_CC, NCOMP_TOTAL, _TOTAL, OPT__REF_FLU_INT_SCHEME,
                                  PhaseUnwrapping_No, Monotonicity,
                                  INT_OPP_SIGN_0TH_ORDER, false );

#  endif // #if ( MODEL == ELBDM ) ... else


// 3.2.2-2. potential
#  ifdef GRAVITY
   const int CSize_Pot3[3] = { CSize_Pot, CSize_Pot, CSize_Pot };
   const int CStart_Pot[3] = { CGhost_Pot, CGhost_Pot, CGhost_Pot };

   real *const CData_Pot = CData_Next;
   CData_Next += CUBE( CSize_Pot );

   real (*FData_Pot)[FSize_CC][FSize_CC] = new real [FSize_CC][FSize_CC][FSize_CC];

   Interpolate( CData_Pot, CSize_Pot3, CStart_Pot, CRange_CC, &FData_Pot[0][0][0],
                FSize_CC3, FStart_CC, 1, OPT__REF_POT_INT_SCHEME, PhaseUnwrapping_No, &Monotonicity_No,
                INT_MONO_COEFF, IntOppSign0thOrder_No );
#  endif


// 3.2.2-3. magnetic field
#  ifdef MHD
   const int FSize_Mag [3][3] = {  { PS2P1, PS2,   PS2   },
                                   { PS2,   PS2P1, PS2   },
                                   { PS2,   PS2,   PS2P1 }  };
   const int FStart_Mag[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
   const int CSize_Mag [3][3] = {  { CSize_Mag_N, CSize_Mag_T, CSize_Mag_T },
                                   { CSize_Mag_T, CSize_Mag_N, CSize_Mag_T },
                                   { CSize_Mag_T, CSize_Mag_T, CSize_Mag_N }  };
   const int CStart_Mag[3][3] = {  { 0,          CGhost_Mag, CGhost_Mag },
                                   { CGhost_Mag, 0,          CGhost_Mag },
                                   { CGhost_Mag, CGhost_Mag, 0          }  };
   const int CRange_Mag[3]    = { PS1, PS1, PS1 };

   real *const CData_MagX = CData_Next + MAGX*CSize_Mag_N*SQR( CSize_Mag_T );
   real *const CData_MagY = CData_Next + MAGY*CSize_Mag_N*SQR( CSize_Mag_T );
   real *const CData_MagZ = CData_Next + MAGZ*CSize_Mag_N*SQR( CSize_Mag_T );
   CData_Next += NCOMP_MAG*CSize_Mag_N*SQR( CSize_Mag_T );

   real (*FData_Mag)[ PS2P1*SQR(PS2) ] = new real [NCOMP_MAG][ PS2P1*SQR(PS2) ];

   const real *CData_Mag3v[NCOMP_MAG] = { CData_MagX, CData_MagY, CData_MagZ };
         real *FData_Mag3v[NCOMP_MAG] = { FData_Mag[MAGX], FData_Mag[MAGY], FData_Mag[MAGZ] };

// set the B field on the coarse-fine interfaces
   const real *Mag_FInterface_Ptr[6] = { NULL, NULL, NULL, NULL, NULL, NULL };

   for (int s=0; s<6; s++)
   {
      const int TRank = CFB_SibRank[s];

//    we set TRank>=0 on the coarse-fine interfaces
      if ( TRank >= 0 )
      {
         Mag_FInterface_Ptr[s] = CFB_BFieldEachRank[TRank] + CFB_OffsetEachRank[TRank];

         CFB_OffsetEachRank[TRank] += SQR( PS2 );
      }
   }

// perform divergence-free interpolation
   MHD_InterpolateBField( CData_Mag3v, CSize_Mag, CStart_Mag, CRange_Mag,
                          FData_Mag3v, FSize_Mag, FStart_Mag, Mag_FInterface_Ptr,
                          OPT__REF_MAG_INT_SCHEME, Monotonicity_Yes );
#  endif // #ifdef MHD


// 3.2.3 check minimum density and pressure/internal energy
// --> note that it's unnecessary to check negative passive scalars thanks to the monotonic interpolation
// --> but we do renormalize passive scalars here
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM  ||  (defined DENS && NCOMP_PASSIVE>0) )
   for (int k=0; k<FSize_CC; k++)
   for (int j=0; j<FSize_CC; j++)
   for (int i=0; i<FSize_CC; i++)
   {
//    check minimum density
      const real DensOld = FData_Flu[DENS][k][j][i];

      if ( DensOld < MIN_DENS )
      {
//       rescale wave function (unnecessary if OPT__INT_PHASE if off, in which case we will rescale all wave functions later)
#        if ( MODEL == ELBDM )
         if ( OPT__INT_PHASE )
         {
            const real Rescale = SQRT( (real)MIN_DENS / DensOld );

            FData_Flu[REAL][k][j][i] *= Rescale;
            FData_Flu[IMAG][k][j][i] *= Rescale;
         }
#        endif

//       apply density floor
         FData_Flu[DENS][k][j][i] = MIN_DENS;
      }


#     if ( MODEL == HYDRO )
//    compute magnetic energy
#     ifdef MHD
      const real Emag = MHD_GetCellCenteredBEnergy( FData_Mag[MAGX], FData_Mag[MAGY], FData_Mag[MAGZ],
                                                    PS2, PS2, PS2, i, j, k );
#     else
      const real Emag = NULL_REAL;
#     endif

//    ensure consistency between pressure, total energy density, and the dual-energy variable
//    --> here we ALWAYS use the dual-energy variable to correct the total energy density
//    --> we achieve that by setting the dual-energy switch to an extremely larger number and ignore
//        the runtime parameter DUAL_ENERGY_SWITCH here
#     ifdef DUAL_ENERGY
      const bool CheckMinPres_Yes = true;
      const real UseEnpy2FixEngy  = HUGE_NUMBER;
      char dummy;    // we do not record the dual-energy status here

      Hydro_DualEnergyFix( FData_Flu[DENS][k][j][i], FData_Flu[MOMX][k][j][i], FData_Flu[MOMY][k][j][i],
                           FData_Flu[MOMZ][k][j][i], FData_Flu[ENGY][k][j][i], FData_Flu[ENPY][k][j][i],
                           dummy, EoS_AuxArray_Flt[1], EoS_AuxArray_Flt[2], CheckMinPres_Yes, MIN_PRES,
                           UseEnpy2FixEngy, Emag );

#     else // #ifdef DUAL_ENERGY

//    apply internal energy floor
      FData_Flu[ENGY][k][j][i]
         = Hydro_CheckMinEintInEngy( FData_Flu[DENS][k][j][i], FData_Flu[MOMX][k][j][i], FData_Flu[MOMY][k][j][i],
                                     FData_Flu[MOMZ][k][j][i], FData_Flu[ENGY][k][j][i], MIN_EINT, Emag );
#     endif // #ifdef DUAL_ENERGY ... else ...
#     endif // #if ( MODEL == HYDRO )


//    normalize passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      if ( OPT__NORMALIZE_PASSIVE )
      {
         real Passive[NCOMP_PASSIVE];

         for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = FData_Flu[ NCOMP_FLUID + v ][k][j][i];

         Hydro_NormalizePassive( FData_Flu[DENS][k][j][i], Passive, PassiveNorm_NVar, PassiveNorm_VarIdx );

         for (int v=0; v<NCOMP_PASSIVE; v++)    FData_Flu[ NCOMP_FLUID + v ][k][j][i] = Passive[v];
      }
#     endif

   } // i,j,k
#  endif // #if ( MODEL == HYDRO  ||  MODEL == ELBDM )


// 3.3 copy data from FData_XXX to patch pointers
   int offset_in[3], i_in, j_in, k_in;

   for (int LocalID=0; LocalID<8; LocalID++)
   {
      const int SonPID = SonPID0 + LocalID;

      offset_in[0] = TABLE_02( LocalID, 'x', 0, PS1 );
      offset_in[1] = TABLE_02( LocalID, 'y', 0, PS1 );
      offset_in[2] = TABLE_02( LocalID, 'z', 0, PS1 );

//    fluid data
      for (int v=0; v<NCOMP_TOTAL; v++)  {
      for (int k=0; k<PS1; k++)  {  k_in = k + offset_in[2];
      for (int j=0; j<PS1; j++)  {  j_in = j + offset_in[1];
      for (int i=0; i<PS1; i++)  {  i_in = i + offset_in[0];

         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[v][k][j][i] = FData_Flu[v][k_in][j_in][i_in];

      }}}}

//    potential data
#     ifdef GRAVITY
      for (int k=0; k<PS1; k++)  {  k_in = k + offset_in[2];
      for (int j=0; j<PS1; j++)  {  j_in = j + offset_in[1];
      for (int i=0; i<PS1; i++)  {  i_in = i + offset_in[0];

         amr->patch[FSg_Pot][SonLv][SonPID]->pot[k][j][i] = FData_Pot[k_in][j_in][i_in];

      }}}
#     endif

//    magnetic field
#     ifdef MHD
      const int Bwidth[3][3] = { {PS1P1, PS1, PS1}, {PS1, PS1P1, PS1}, {PS1, PS1, PS1P1} };
      int idx_B_in, idx_B_out;

      for (int v=0; v<NCOMP_MAG; v++)     {  idx_B_out = 0;
      for (int k=0; k<Bwidth[v][2]; k++)  {  k_in      = k + offset_in[2];
      for (int j=0; j<Bwidth[v][1]; j++)  {  j_in      = j + offset_in[1];
                                             idx_B_in  = IDX321( offset_in[0], j_in, k_in, FSize_Mag[v][0], FSize_Mag[v][1] );
      for (int i=0; i<Bwidth[v][0]; i++)  {

         amr->patch[FSg_Mag][SonLv][SonPID]->magnetic[v][ idx_B_out ++ ] = FData_Mag[v][ idx_B_in ++ ];

      }}}}
#     endif

//    rescale real and imaginary parts to get the correct density in ELBDM if OPT__INT_PHASE is off
#     if ( MODEL == ELBDM )
      real Real, Imag, Rho_Corr, Rho_Wrong, Rescale;

      if ( !OPT__INT_PHASE )
      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {
         Real      = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[REAL][k][j][i];
         Imag      = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = amr->patch[FSg_Flu][SonLv][SonPID]->fluid[DENS][k][j][i];

//       be careful about the negative density introduced from the round-off errors
         if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
         {
            amr->patch[FSg_Flu][SonLv][SonPID]->fluid[DENS][k][j][i] = (real)0.0;
            Rescale = (real)0.0;
         }
         else
            Rescale = SQRT( Rho_Corr/Rho_Wrong );

         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[REAL][k][j][i] *= Rescale;
         amr->patch[FSg_Flu][SonLv][SonPID]->fluid[IMAG][k][j][i] *= Rescale;
      }
#     endif
   } // for (int LocalID=0; LocalID<8; LocalID++)


// 4. pass particles from father to son if they are in the same rank
//    --> otherwise these particles will be transferred to the real son patches by calling
//        Par_PassParticle2Son_MultiPatch() in LB_Refine()
#  ifdef PARTICLE
   if ( FaPID >= 0  &&  FaPID < amr->NPatchComma[FaLv][1] )    Par_PassParticle2Son_SinglePatch( FaLv, FaPID );
#  endif


// free memory
   if ( FaIsHome )   delete [] CData;
   delete [] FData_Flu;
#  ifdef GRAVITY
   delete [] FData_Pot;
#  endif
#  ifdef MHD
   delete [] FData_Mag;
#  endif


   return SonPID0;

} // FUNCTION : AllocateSonPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  DeallocateSonPatch
// Description :  Deallocate eight son patches at FaLv+1
//
// Note        :  1. Just to avoid duplicate code segment
//                2. List NewSonPID0_Real[] will be reset
//
// Parameter   :  FaLv            : Target refinement level to be refined
//                FaPID           : Father patch index to remove sons
//                NNew_Real0      : Number of newly-allocated real patches with LocalID==0
//                NewSonPID0_Real : List recording the indices of all newly-allocated real patches
//                                  with LocalID==0
//                SwitchIdx       : Element in NewSonPID0_Real to be reset
//
//                PARTICLE-only parameters (call-by-reference)
//                RefineS2F_Send_NPatchTotal : Total number of patches for exchanging particles from sons to fathers
//                RefineS2F_Send_PIDList     : Patch indices for exchanging particles from sons to fathers
//-------------------------------------------------------------------------------------------------------
void DeallocateSonPatch( const int FaLv, const int FaPID, const int NNew_Real0, int NewSonPID0_Real[],
                         int SwitchIdx, int &RefineS2F_Send_NPatchTotal, int *&RefineS2F_Send_PIDList )
{

   const int MirSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   const int SonLv      = FaLv + 1;
   const int GraLv      = FaLv + 2;
   const int SonPID0    = amr->patch[0][FaLv][FaPID]->son;
   const int NewPID0    = SonPID0;
   const int OldPID0    = amr->num[SonLv] - 8;

   int NewPID, OldPID, OldGraPID0, OldFaPID, SibPID, MaxSonPID0, MaxIdx, TempSonPID0;


// check : son patch must exist and be real patches
#  ifdef GAMER_DEBUG
   if ( SonPID0 < 0  ||  SonPID0 >= amr->NPatchComma[SonLv][1] )
      Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, SonNReal %d, incorrect SonPID0 = %d !!\n",
                 FaLv, FaPID, amr->NPatchComma[SonLv][1], SonPID0 );
#  endif


// pass particles from sons to father
#  ifdef PARTICLE
   Par_PassParticle2Father( FaLv, FaPID );

// record the father patch indices for exchanging particles if they are buffer patches with particles
// --> these particles will be transferred to the real father patches by calling
//     Par_LB_Refine_SendParticle2Father() in LB_Refine()
   if ( FaPID >= amr->NPatchComma[FaLv][1]  &&  amr->patch[0][FaLv][FaPID]->NPar > 0 )
   {
#     ifdef DEBUG_PARTICLE
      const int FaNBuff = amr->NPatchComma[FaLv][3] - amr->NPatchComma[FaLv][1];
      if ( RefineS2F_Send_NPatchTotal >= FaNBuff )
         Aux_Error( ERROR_INFO, "target index (%d) >= FaNBuff (%d) !!\n", RefineS2F_Send_NPatchTotal, FaNBuff );
#     endif

      RefineS2F_Send_PIDList[ RefineS2F_Send_NPatchTotal ++ ] = FaPID;
   }
#  endif // #ifdef PARTICLE


// deallocate the unflagged child patches and reset sibling indices to -1
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
   {
      for (int s=0; s<26; s++)
      {
         SibPID = amr->patch[0][SonLv][SonPID]->sibling[s];
         if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = -1;
      }

      amr->pdelete( SonLv, SonPID, OPT__REUSE_MEMORY );
   }

// record NPatchComma
   amr->NPatchComma[SonLv][1] -= 8;

// construct relation : father -> son
   amr->patch[0][FaLv][FaPID]->son = -1;

// relink the child patch pointers so that no patch indices are skipped
   if ( NewPID0 != OldPID0 )
   {
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         NewPID = NewPID0 + LocalID;
         OldPID = OldPID0 + LocalID;

//       swap patch pointers between the old and new PID
//       --> works no matter OPT__REUSE_MEMORY is on or off
//       --> note that when OPT__REUSE_MEMORY is on, we don't want to set pointers of OldPID as NULL
         for (int Sg=0; Sg<2; Sg++)
            Aux_SwapPointer( (void**)&amr->patch[Sg][SonLv][OldPID], (void**)&amr->patch[Sg][SonLv][NewPID] );

//       swap back fluid[] with 1-FluSg, pot[] with 1-PotSg, and magnetic[] with 1-MagSg since
//       we use them to temporarily store the buffer patch data if OPT__REUSE_MEMORY is used
//       --> see the description in LB_Refine_AllocateNewPatch
//       --> ugly trick ...
         if ( OPT__REUSE_MEMORY )
         {
            const int FSg_Flu2 = 1 - amr->FluSg[SonLv];
            Aux_SwapPointer( (void**)&amr->patch[FSg_Flu2][SonLv][OldPID]->fluid,
                             (void**)&amr->patch[FSg_Flu2][SonLv][NewPID]->fluid );

#           ifdef GRAVITY
            const int FSg_Pot2 = 1 - amr->PotSg[SonLv];
            Aux_SwapPointer( (void**)&amr->patch[FSg_Pot2][SonLv][OldPID]->pot,
                             (void**)&amr->patch[FSg_Pot2][SonLv][NewPID]->pot );
#           endif

#           ifdef MHD
            const int FSg_Mag2 = 1 - amr->MagSg[SonLv];
            Aux_SwapPointer( (void**)&amr->patch[FSg_Mag2][SonLv][OldPID]->magnetic,
                             (void**)&amr->patch[FSg_Mag2][SonLv][NewPID]->magnetic );
#           endif
         }

//       reconstruct relation : grandson -> son
         OldGraPID0 = amr->patch[0][SonLv][NewPID]->son;
         if ( OldGraPID0 >= 0 )
         {
#           ifdef GAMER_DEBUG
            if ( GraLv >= NLEVEL )
               Aux_Error( ERROR_INFO, "SonLv %d, NewPID %d, OldGraPID0 %d, GraLv %d >= NLEVEL (%d) !!\n",
                          SonLv, NewPID, OldGraPID0, GraLv, NLEVEL );
#           endif

            for (int OldGraPID=OldGraPID0; OldGraPID<OldGraPID0+8; OldGraPID++)
               amr->patch[0][GraLv][OldGraPID]->father = NewPID;
         }

//       reconstruct relation : sibling
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][SonLv][NewPID]->sibling[s];
            if ( SibPID >= 0 )   amr->patch[0][SonLv][SibPID]->sibling[ MirSib[s] ] = NewPID;
         }
      } // for (int LocalID=0; LocalID<8; LocalID++)

//    reconstruct relation : father -> son
      OldFaPID = amr->patch[0][SonLv][NewPID0]->father;
#     ifdef GAMER_DEBUG
      if ( OldFaPID == -1 )
         Aux_Error( ERROR_INFO, "SonLv %d, NewPID0 %d -> father == -1 !!\n", SonLv, NewPID0 );
#     endif
      amr->patch[0][FaLv][OldFaPID]->son = NewPID0;

//    reset the list "NewSonPID0_Real"
      if ( NNew_Real0 > 0 )
      {
         if ( SwitchIdx >= 0 )   // not all elements in NewSonPID0_Real have been reset
         {
#           ifdef GAMER_DEBUG
            if ( NewSonPID0_Real[SwitchIdx] != OldPID0 )
               Aux_Error( ERROR_INFO, "SonLv %d, NewSonPID0_Real[%d] (%d) != OldPID0 (%d) !!\n",
                          SonLv, SwitchIdx, NewSonPID0_Real[SwitchIdx], OldPID0 );
#           endif

            NewSonPID0_Real[SwitchIdx] = NewPID0;

//          put the largest element on "NewSonPID0_Real[NNew_Real0-1]"
            if ( SwitchIdx-1 < 0 )
            {
               MaxSonPID0 = -__INT_MAX__;
               MaxIdx     = -__INT_MAX__;

               for (int t=0; t<NNew_Real0; t++)
               {
                  if ( NewSonPID0_Real[t] > MaxSonPID0 )
                  {
                     MaxSonPID0 = NewSonPID0_Real[t];
                     MaxIdx     = t;
                  }
               }

               TempSonPID0                   = NewSonPID0_Real[NNew_Real0-1];
               NewSonPID0_Real[NNew_Real0-1] = MaxSonPID0;
               NewSonPID0_Real[MaxIdx      ] = TempSonPID0;
            }
         } // if ( SwitchIdx >= 0 )

         else  // all elements in NewSonPID0_Real have been reset --> need to resort the list
         {
            SwitchIdx = NNew_Real0 - 1;

            if ( NewSonPID0_Real[SwitchIdx] == OldPID0 )
            {
               NewSonPID0_Real[SwitchIdx] = NewPID0;

//             put the largest element on "NewSonPID0_Real[NNew_Real0-1]"
               MaxSonPID0 = -__INT_MAX__;
               MaxIdx     = -__INT_MAX__;

               for (int t=0; t<NNew_Real0; t++)
               {
                  if ( NewSonPID0_Real[t] > MaxSonPID0 )
                  {
                     MaxSonPID0 = NewSonPID0_Real[t];
                     MaxIdx     = t;
                  }
               }

               TempSonPID0                   = NewSonPID0_Real[NNew_Real0-1];
               NewSonPID0_Real[NNew_Real0-1] = MaxSonPID0;
               NewSonPID0_Real[MaxIdx      ] = TempSonPID0;
            }
         } // if ( SwitchIdx >= 0 ) ... else ...
      } // if ( NNew_Real0 > 0 )
   } // if ( NewPID0 != OldPID0 )

} // FUNCTION : DeallocateSonPatch



#endif // #ifdef LOAD_BALANCE
