#include "GAMER.h"
#include "ReadPara.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_LoadExtPotTable
// Description :  Load the external potential table
//
// Note        :  1. Invoked by Init_ExtAccPot()
//                2. Enabled by the runtime option "OPT__EXT_POT=EXT_POT_TABLE"
//                3. Floating-point type is set by FLOAT8 for now
//                   --> EXT_POT_TABLE_FLOAT8 is NOT supported yet
//                4. The loaded table will be sent to GPU by invoking CUAPI_SendExtPotTable2GPU()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_LoadExtPotTable()
{

// only for EXT_POT_TABLE
   if ( OPT__EXT_POT != EXT_POT_TABLE )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// checks
// parameters
   for (int d=0; d<3; d++) {
      if ( EXT_POT_TABLE_NPOINT[d] < 2 )
         Aux_Error( ERROR_INFO, "EXT_POT_TABLE_NPOINT[%d] = %d < 2 !!\n", d, EXT_POT_TABLE_NPOINT[d] );
   }

   if ( EXT_POT_TABLE_DH <= 0.0 )
      Aux_Error( ERROR_INFO, "EXT_POT_TABLE_DH = %14.7e <= 0.0 !!\n", EXT_POT_TABLE_DH );

#  ifdef FLOAT8
   if ( ! EXT_POT_TABLE_FLOAT8 )
      Aux_Error( ERROR_INFO, "must enable EXT_POT_TABLE_FLOAT8 since FLOAT8 is on !!\n" );
#  else
   if (   EXT_POT_TABLE_FLOAT8 )
      Aux_Error( ERROR_INFO, "must disable EXT_POT_TABLE_FLOAT8 since FLOAT8 is off !!\n" );
#  endif

// table must cover the simulation domain plus GRA_GHOST_SIZE (GRA_GHOST_SIZE-0.5, to be more specific)
// level-one cells on each side
// --> for CPU_ExtPotSolver() to fill in the ghost-zone potential (even when STORE_POT_GHOST is off)
//
//  |: Cell boundary
//  ‖: Box boundary
//  ⊚: Cell-centered value in non-ghost-zone
//  ∆: Cell-centered value in ghost-zone
//
//  EXT_POT_TABLE_NPOINT : Total number of ⊚ and ∆
//
//------------------------------------------------------------
//  EXT_POT_TABLE_EDGEL                 EXT_POT_TABLE_EDGER 
//         ↓                                   ↓
//       |-∆-|-∆-‖-⊚-|-⊚-|-⊚-|-⊚-|-⊚-|-⊚-‖-∆-|-∆-|
//               ↑                       ↑ 
//          amr->BoxEdgeL           amr->BoxEdgeR 
//------------------------------------------------------------
//
   const double ExtraWidth = 0.5*GRA_GHOST_SIZE*amr->dh[0]; // use GRA_GHOST_SIZE instead of GRA_GHOST_SIZE-0.5 for extra safety

   for (int d=0; d<3; d++) {
      if ( EXT_POT_TABLE_EDGEL[d] == NoDef_double )
         Aux_Error( ERROR_INFO, "EXT_POT_TABLE_EDGEL[%d] has not been set !!\n", d );

      if ( EXT_POT_TABLE_EDGEL[d] > amr->BoxEdgeL[d]-ExtraWidth )
         Aux_Error( ERROR_INFO, "incorrect left boundary of the external potential table (table=%13.7e > min=%13.7e) !!\n",
                    EXT_POT_TABLE_EDGEL[d], amr->BoxEdgeL[d]-ExtraWidth );

      if ( EXT_POT_TABLE_EDGEL[d]+(EXT_POT_TABLE_NPOINT[d]-1)*EXT_POT_TABLE_DH < amr->BoxEdgeR[d]+ExtraWidth )
         Aux_Error( ERROR_INFO, "incorrect right boundary of the external potential table (table=%13.7e < max=%13.7e) !!\n",
                    EXT_POT_TABLE_EDGEL[d]+(EXT_POT_TABLE_NPOINT[d]-1)*EXT_POT_TABLE_DH, amr->BoxEdgeR[d]+ExtraWidth );
   }

// file existence
   if ( ! Aux_CheckFileExist(EXT_POT_TABLE_NAME) )
      Aux_Error( ERROR_INFO, "external potential table \"%s\" does not exist !!\n", EXT_POT_TABLE_NAME );

// file size
   FILE *FileTemp = fopen( EXT_POT_TABLE_NAME, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long NPoint3D   = (long)EXT_POT_TABLE_NPOINT[0]*EXT_POT_TABLE_NPOINT[1]*EXT_POT_TABLE_NPOINT[2];
   const long ExpectSize = NPoint3D*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the external potential table <%s> (%ld) != expect (%ld) !!\n",
                 EXT_POT_TABLE_NAME, FileSize, ExpectSize );

   fclose( FileTemp );

// memory allocation
   if ( h_ExtPotTable == NULL )
      Aux_Error( ERROR_INFO, "h_ExtPotTable[] has not been allocated !!\n" );

   MPI_Barrier( MPI_COMM_WORLD );


// load table to CPU
   FILE *File = fopen( EXT_POT_TABLE_NAME, "rb" );
   fread( h_ExtPotTable, sizeof(real), NPoint3D, File );
   fclose( File );

// normalized external potential by UNIT_D*UNIT_V**2
   if ( OPT__UNIT )
      for (long int Idx=0; Idx<NPoint3D; Idx++  )
          h_ExtPotTable[Idx] /= UNIT_D * UNIT_V * UNIT_V;

// transfer table to GPU
#  ifdef GPU
   CUAPI_SendExtPotTable2GPU( h_ExtPotTable );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_LoadExtPotTable



#endif // #ifdef GRAVITY
