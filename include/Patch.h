#ifndef __PATCH_H__
#define __PATCH_H__



#include "Macro.h"

#ifdef PARTICLE
#  include <math.h>
#  include "Particle.h"
#endif


void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
ulong Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
long  LB_Corner2Index( const int lv, const int Corner[], const Check_t Check );




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t
// Description :  Data structure of a single patch
//
// Data Member :  fluid           : Fluid variables (mass density, momentum density x, y ,z, energy density)
//                                  --> Including passively advected variables (e.g., metal density)
//                magnetic        : Magnetic field (Bx, By, Bz)
//                pot             : Potential
//                pot_ext         : Potential with GRA_GHOST_SIZE ghost cells on each side
//                                  --> Allocated only if STORE_POT_GHOST is on
//                                  --> Ghost-zone potential are obtained from the Poisson solver directly
//                                      (not from exchanging potential between sibling patches)
//                                  --> Currently it is used for Par->ImproveAcc and SF_CreateStar_AGORA() only
//                                  --> Currently it's useless for buffer patches
//                de_status       : Assigned to (DE_UPDATED_BY_ETOT / DE_UPDATED_BY_DUAL / DE_UPDATED_BY_MIN_PRES /
//                                               DE_UPDATED_BY_ETOT_GRA)
//                                  to indicate whether each cell is updated by the total energy, dual energy variable,
//                                  or minimum allowed pressure
//                                  --> DE_UPDATED_BY_XXX are defined in Macro.h
//                                  --> It's a character array with the size PS1^3
//                                  --> Currently it's only allocated for Sg=0
//                rho_ext         : Density with RHOEXT_GHOST_SIZE (typically 2) ghost cells on each side
//                                  --> Only allocated **temporarily** in the function Prepare_PatchData for storing
//                                      particle mass density
//                                      --> Also note that this array does NOT necessary store the correct particle mass density
//                                          (especially for cells adjacent to the C-C and C-F boundaries) and thus should
//                                          NOT be used outside Prepare_PatchData)
//                                  --> Note that even with NGP mass assignment which requires no ghost zone we
//                                      still allocate rho_ext as (PS1+RHOEXT_GHOST_SIZE)^3
//                flux[6]         : Fluid flux (for the flux-correction operation)
//                                  --> Including passively advected flux (for the flux-correction operation)
//                flux_tmp[6]     : Temporary fluid flux for the option "AUTO_REDUCE_DT"
//                flux_bitrep[6]  : Fluid flux for achieving bitwise reproducibility (i.e., ensuring that the round-off errors are
//                                  exactly the same in different parallelization parameters/strategies)
//                electric        : Electric field for the MHD fix-up operation
///                                 --> Array structure = [sibling index][E field index][cell index]
//                                  --> For sibling indices 0-5, there are two E fields on each face
//                                      --> E field index on x faces: [0/1] = Ey/Ez
//                                                           y faces: [0/1] = Ez/Ex
//                                                           z faces: [0/1] = Ex/Ey
//                                          Cell index dimension on x faces: [Nz][Ny] (= [PS1-1][PS1] for Ey and [PS1][PS1-1] for Ez)
//                                                     dimension on y faces: [Nx][Nz] (= [PS1-1][PS1] for Ez and [PS1][PS1-1] for Ex)
//                                                     dimension on z faces: [Ny][Nx] (= [PS1-1][PS1] for Ex and [PS1][PS1-1] for Ey)
//                                  --> For sibling indices 6-17, there is only one E field on each edge
//                                      --> E field index is always 0 (6-9->Ez, 10-13->Ex, 14-17->Ey)
//                                          Cell index dimension is always [PS1]
//                                  --> To unify the array structure of sibling indices 0-5 and 6-17, we actually convert the
//                                      2D array [E field index][cell index] into a 1D array
//                electric_tmp    : Temporary electric field for the option "AUTO_REDUCE_DT"
//                electric_bitrep : Electric field for achieving bitwise reproducibility in MHD (i.e., ensuring that the
//                                  round-off errors are exactly the same in different parallelization parameters/strategies)
//                ele_corrected   : Array recording whether each component in electric[] has been corrected by the fix-up operation
//                corner[3]       : Grid indices of the cell at patch corner
//                                  --> Note that for an external patch its recorded "corner" will lie outside
//                                      the simulation domain. In other words, periodicity is NOT used to
//                                      convert corner to that of the corresponding real patch
//                                  --> For example, external patches can have corner < 0
//                                  --> Different from EdgeL/R, which always assume periodicity
//                sibling[26]     : Patch IDs of the 26 sibling patches (-1->no sibling; -1XX->external)
//
//                                  NOTE FOR NON-PERIODIC BOUNDARY CONDITIONS:
//                                  If a target sibling patch is an external patch (which lies outside the simulation domain),
//                                  the corresponding sibling index is set to "SIB_OFFSET_NONPERIODIC-Sibling", where Sibling
//                                  represents the sibling direction of the boundary region.
//                                  (e.g., if amr->sibling[XX] lies outside the -y boundary (boundary sibling = 2), we set
//                                   amr->sibling[XX] = SIB_OFFSET_NONPERIODIC-2 = -102 for SIB_OFFSET_NONPERIODIC==-100)
//                                  --> All patches without sibling patches along the XX direction will have
//                                      "amr->sibling[XX] < 0"
//
//                father          : Patch ID of the father patch
//                son             : Patch ID of the child patch (-1->no son)
//
//                                  LOAD_BALANCE NOTE:
//                                  For patches with sons living abroad (not at the same rank as iteself), their son indices
//                                  are set to "SON_OFFSET_LB-SonRank", where SonRank represents the MPI rank where their
//                                  sons live. (e.g., if the real son at rank 123, the son index is set to
//                                  SON_OFFSET_LB-123=-1123 for SON_OFFSET_LB==-1000)
//                                  --> All patches (i.e., both real and buffer patches) with sons will have SonPID != -1
//
//                flag            : Refinement flag (true/false)
//                Active          : Used by OPT__REUSE_MEMORY to indicate whether this patch is active or inactive
//                                  --> active:    patch has been allocated and activated   (included in   num[lv])
//                                      inactive:  patch has been allocated but deactivated (excluded from num[lv])
//                                  --> Note that active/inactive have nothing to do with the allocation of field arrays (e.g., fluid)
//                                      --> For both active and inactive patches, field arrays may be allocated or == NULL
//                                  --> However, currently the flux arrays (i.e., flux, flux_tmp, and flux_bitrep) are guaranteed
//                                      to be NULL for inactive patches
//                EdgeL/R         : Left and right edge of the patch
//                                  --> Note that we always apply periodicity to EdgeL/R. So for an external patch its
//                                      recorded "EdgeL/R" will still lie inside the simulation domain and will be
//                                      exactly the same as the EdgeL/R of the corresponding real patches.
//                                  --> For example, the external patches just outside the simulation left edge will have
//                                      EdgeL = BoxEdgeR-PatchSize*dh[lv] and EdgeR = BoxEdgeR, and for those just outside
//                                      the simulation right edge will have EdgeL = BoxEdgeL and EdgeR = BoxEdgeL+PatchSize*dh[lv]
//                                  --> Different from corner[3], which do NOT assume periodicity
//                PaddedCr1D      : 1D corner coordiniate padded with two base-level patches on each side
//                                  in each direction, normalized to the finest-level patch scale (PATCH_SIZE)
//                                  --> Each PaddedCr1D defines a unique 3D position
//                                  --> Patches at different levels with the same PaddedCr1D have the same
//                                      3D corner coordinates
//                                  --> This number is independent of periodicity (because of the padded patches)
//                LB_Idx          : Space-filling-curve index for load balance
//                NPar            : Number of particles belonging to this leaf patch
//                ParListSize     : Size of the array ParList (ParListSize can be >= NPar)
//                ParList         : List recording the IDs of all particles belonging to this leaf real patch
//                NPar_Copy       : Number of particles collected from other patches. There are three kinds of patches that
//                                  can have NPar_Copy != -1
//                                  --> (1) Non-leaf real   patches (both SERIAL and LOAD_BALANCE)
//                                          --> Particles collected from the descendants (sons, grandsons, ...) of this patch
//                                      (2) Non-leaf buffer patches (LOAD_BALANCE only)
//                                          --> Particles collected from the corresponding non-leaf real patches in (1)
//                                      (3) Leaf     buffer patches (LOAD_BALANCE only)
//                                          --> Particles collected from the corresponding leaf real patches
//                                              (those with NPar and ParList property set)
//                                  --> **Leaf real** patches will always have NPar_Copy == -1
//                                  --> In SERIAL mode, these non-leaf real patches will have their particle
//                                      IDs stored in ParList_Copy, which points to the same particle repository
//                                      (i.e., the amr->Par->Attribute[]). In comparison, in LOAD_BALANCE mode,
//                                      since particles corresponding to NPar_Copy may be collected from other ranks,
//                                      these patches will allocate a local particle attribute array called ParMassPos_Copy
//                                      (since currently we only collect mass and position for these particles).
//                                  --> Note that non-leaf patches may have NPar>0 temporarily after updating particle position.
//                                      It's because particles travelling from coarse to fine grids will stay in coarse grids
//                                      temporarily until the velocity correction is done.
//                                      --> For these patches, NPar_Copy will be **the sum of NPar and the number of particles
//                                          collected from other patches**, and ParList_Copy (or ParMassPos_Copy) will contain
//                                          information of particles belonging to NPar as well.
//                                      --> It makes implementation simplier. For leaf real patches, one only needs to consider
//                                          NPar and ParList. While for all other patches, one only needs to consider NPar_Copy,
//                                          ParList_Copy (or ParMassPos_Copy). One never needs to consider both.
//                ParList_Copy    : List recording the IDs of all particles belonging to the descendants of this patch
//                                  (and particles temporarily locate in this patch waiting for the velocity correction, see
//                                  discussion above)
//                                  --> for SERIAL only
//                ParMassPos_Copy : Pointer arrays storing the mass and position of NPar_Copy particles collected from other patches
//                                  --> for LOAD_BALANCE only
//                NPar_Escp       : Number of particles escaping from this patch
//                ParList_Escp    : List recording the IDs of all particles escaping from this patch
//
// Method      :  patch_t         : Constructor
//               ~patch_t         : Destructor
//                Activate        : Activate patch
//                fnew            : Allocate flux[]
//                fdelete         : Deallocate flux[]
//                enew            : Allocate electric[]
//                edelete         : Deallocate electric[]
//                hnew            : Allocate fluid[]
//                hdelete         : Deallocate fluid[]
//                mnew            : Allocate magnetic[]
//                mdelete         : Deallocate magnetic[]
//                gnew            : Allocate pot[]
//                gdelete         : Deallocate pot[]
//                snew            : Allocate de_status[]
//                sdelete         : Deallocate de_status[]
//                dnew            : Allocate rho_ext[]
//                ddelete         : Deallocate rho_ext[]
//                AddParticle     : Add particles to the particle list
//                RemoveParticle  : Remove particles from the particle list
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid)[PS1][PS1][PS1];

#  ifdef MHD
   real (*magnetic)[ PS1P1*SQR(PS1) ];
#  endif

#  ifdef GRAVITY
   real (*pot)[PS1][PS1];
#  ifdef STORE_POT_GHOST
   real (*pot_ext)[GRA_NXT][GRA_NXT];
#  endif
#  endif // GRAVITY

#  ifdef DUAL_ENERGY
   char (*de_status)[PS1][PS1];
#  endif

#  ifdef PARTICLE
   real (*rho_ext)[RHOEXT_NXT][RHOEXT_NXT];
#  endif

   real (*flux       [6])[PS1][PS1];
   real (*flux_tmp   [6])[PS1][PS1];
#  ifdef BIT_REP_FLUX
   real (*flux_bitrep[6])[PS1][PS1];
#  endif

#  ifdef MHD
   real (*electric       [18]);
   real (*electric_tmp   [18]);
#  ifdef BIT_REP_ELECTRIC
   real (*electric_bitrep[18]);
#  endif
   bool ele_corrected[12];
#  endif

   int    corner[3];
   int    sibling[26];
   int    father;
   int    son;
   bool   flag;
   bool   Active;
   double EdgeL[3];
   double EdgeR[3];

   ulong  PaddedCr1D;
   long   LB_Idx;

#  ifdef PARTICLE
   int    NPar;
   int    NParType[2];
   int    ParListSize;
   long  *ParList;

   int    NPar_Copy;
#  ifdef LOAD_BALANCE
   real  *ParMassPos_Copy[4];
#  else
   long  *ParList_Copy;
#  endif

   int    NPar_Escp[26];
   long  *ParList_Escp[26];
#  endif



   //===================================================================================
   // Constructor :  patch_t
   // Description :  Constructor of the structure "patch_t"
   //
   // Note        :  Initialize data members
   //
   // Parameter   :  scale_x/y/z : Grid scale indices of the patch corner
   //                FaPID       : Patch ID of the father patch
   //                FluData     : true --> Allocate the hydrodynamic array fluid[]
   //                                       --> rho_ext[] will NOT be allocated here even if PARTICLE is on
   //                MagData     : true --> Allocate the MHD array magnetic[]
   //                                       --> Useless if "MHD" is turned off
   //                PotData     : true --> Allocate the potential array pot[] (has no effect if "GRAVITY" is turned off)
   //                                       --> pot_ext[] will be allocated as well if STORE_POT_GHOST is on
   //                DE_Status   : true --> Allocate the dual-energy status array de_status[]
   //                                       --> Useless if "DUAL_ENERGY" is turned off
   //                lv          : Refinement level of the newly created patch
   //                BoxScale    : Simulation box scale
   //                BoxEdgeL    : Simulation box left edge
   //                dh_min      : Cell size at the maximum level
   //===================================================================================
   patch_t( const int scale_x, const int scale_y, const int scale_z, const int FaPID, const bool FluData,
            const bool MagData, const bool PotData, const bool DE_Status, const int lv, const int BoxScale[],
            const double BoxEdgeL[], const double dh_min )
   {

//    always initialize field pointers (e.g., fluid, pot, ...) as NULL if they are not allocated here
      const bool InitPtrAsNull_Yes = true;
      Activate( scale_x, scale_y, scale_z, FaPID, FluData, MagData, PotData, DE_Status, lv, BoxScale,
                BoxEdgeL, dh_min, InitPtrAsNull_Yes );

   } // METHOD : patch_t



   //===================================================================================
   // Constructor :  Activate
   // Description :  Initialize data members
   //
   // Note        :  Called by the patch constructor "patch_t" and amr->pnew when OPT__REUSE_MEMORY is adopted
   //
   // Parameter   :  scale_x/y/z   : Grid scale indices of the patch corner
   //                FaPID         : Patch ID of the father patch
   //                FluData       : true --> Allocate the hydrodynamic array fluid[]
   //                                         --> rho_ext[] will NOT be allocated here even if PARTICLE is on
   //                MagData       : true --> Allocate the MHD array magnetic[]
   //                                         --> Useless if "MHD" is turned off
   //                PotData       : true --> Allocate the potential array pot[] (has no effect if "GRAVITY" is turned off)
   //                                         --> pot_ext[] will be allocated as well if STORE_POT_GHOST is on
   //                DE_Status     : true --> Allocate the dual-energy status array de_status[]
   //                                         --> Useless if "DUAL_ENERGY" is turned off
   //                lv            : Refinement level of the newly created patch
   //                BoxScale      : Simulation box scale
   //                BoxEdgeL      : Simulation box left edge
   //                dh_min        : Cell size at the maximum level
   //                InitPtrAsNull : Whether or not to initialize the field arrays
   //                                (i.e., fluid, pot, pot_ext, rho_ext, de_status) as NULL
   //                                --> It is used mainly for OPT__REUSE_MEMORY, where we don't want to set these
   //                                    pointers as NULL if the patch has been allocated but marked as inactive
   //                                    --> Since the field arrays may be allocated already and we want to reuse them
   //                                --> But note that we must initialize these pointers as NULL when allocating
   //                                    (not activating) patches and before calling hnew and gnew
   //                                    --> otherwise these pointers become ill-defined, which will make hdelete and
   //                                        gdelete crash
   //                                --> Does NOT apply to flux arrays (i.e., flux, flux_tmp, and flux_bitrep) which
   //                                    are always initialized as NULL here
   //                                --> Does not apply to any particle variable (except rho_ext)
   //===================================================================================
   void Activate( const int scale_x, const int scale_y, const int scale_z, const int FaPID, const bool FluData,
                  const bool MagData, const bool PotData, const bool DE_Status, const int lv, const int BoxScale[],
                  const double BoxEdgeL[], const double dh_min, const bool InitPtrAsNull )
   {

      corner[0] = scale_x;
      corner[1] = scale_y;
      corner[2] = scale_z;
      father    = FaPID;
      son       = -1;
      flag      = false;
      Active    = true;

      for (int s=0; s<26; s++ )  sibling[s] = -1;     // -1 <--> NO sibling

      const int Padded              = 1<<NLEVEL;
      const int BoxNScale_Padded[3] = { BoxScale[0]/PS1 + 2*Padded,
                                        BoxScale[1]/PS1 + 2*Padded,
                                        BoxScale[2]/PS1 + 2*Padded }; // normalized and padded box scale
      int Cr_Padded[3];
      for (int d=0; d<3; d++)    Cr_Padded[d] = corner[d]/PS1 + Padded;

#     ifdef GAMER_DEBUG
      for (int d=0; d<3; d++)
      {
         if ( Cr_Padded[d] < 0  ||  Cr_Padded[d] >= BoxNScale_Padded[d] )
            Aux_Error( ERROR_INFO, "incorrect Cr_Padded[%d]=(%d) --> out of range [%d ... %d) !!\n",
                       d, Cr_Padded[d], 0, BoxNScale_Padded[d] );
      }
#     endif

      PaddedCr1D = Mis_Idx3D2Idx1D( BoxNScale_Padded, Cr_Padded );   // independent of periodicity
      LB_Idx     = LB_Corner2Index( lv, corner, CHECK_OFF );         // always assumes periodicity

//    set the patch edge
      const int PScale = PS1*( 1<<(TOP_LEVEL-lv) );
      for (int d=0; d<3; d++)
      {
//       to ensure that buffer patches and their corresponding real patches have the same EdgeL/R, do not write EdgeR[d] as
//       EdgeR[d] = BoxEdgeL[d] + (double)(  ( corner[d] + BoxScale[d] + PScale ) % BoxScale[d] )*dh_min;
//       --> otherwise the buffer patches just outside the simulation left edge (and the real patches just inside the simulation
//           right edge) will have EdgeR=BoxEdgeL instead of EdgeR=BoxEdgeR
//       --> assuming periodicity
         EdgeL[d] = BoxEdgeL[d] + (double)(  ( corner[d] + BoxScale[d] ) % BoxScale[d]           )*dh_min;
         EdgeR[d] = BoxEdgeL[d] + (double)(  ( corner[d] + BoxScale[d] ) % BoxScale[d] + PScale  )*dh_min;

//       do no use the following non-periodic version anymore --> it does not work with the current particle implementation
         /*
         EdgeL[d] = BoxEdgeL[d] + (double)corner[d]*dh_min;
         EdgeR[d] = BoxEdgeL[d] + (double)( corner[d] + PScale )*dh_min;
         */
      }

//    must initialize these pointers as NULL when allocating (not activating) patches and before calling hnew and gnew
//    --> otherwise these pointers become ill-defined, which will make hdelete and gdelete crash
      if ( InitPtrAsNull )
      {
         fluid     = NULL;

#        ifdef MHD
         magnetic  = NULL;
#        endif

#        ifdef GRAVITY
         pot       = NULL;
#        ifdef STORE_POT_GHOST
         pot_ext   = NULL;
#        endif
#        endif // #ifdef GRAVITY

#        ifdef DUAL_ENERGY
         de_status = NULL;
#        endif

#        ifdef PARTICLE
         rho_ext   = NULL;
#        endif
      }

      for (int s=0; s<6; s++)
      {
         flux       [s] = NULL;
         flux_tmp   [s] = NULL;
#        ifdef BIT_REP_FLUX
         flux_bitrep[s] = NULL;
#        endif
      }

#     ifdef MHD
      for (int s=0; s<18; s++)
      {
         electric       [s] = NULL;
         electric_tmp   [s] = NULL;
#        ifdef BIT_REP_ELECTRIC
         electric_bitrep[s] = NULL;
#        endif
      }
#     endif

      if ( FluData )    hnew();
#     ifdef MHD
      if ( MagData )    mnew();
#     endif
#     ifdef GRAVITY
      if ( PotData )    gnew();
#     endif
#     ifdef DUAL_ENERGY
      if ( DE_Status )  snew();
#     endif

#     ifdef PARTICLE
      NPar         = 0;          // must be initialized as 0
      NParType[0]  = 0;          // must be initialized as 0
      NParType[1]  = 0;          // must be initialized as 0
      ParListSize  = 0;          // must be initialized as 0
      ParList      = NULL;

      NPar_Copy    = -1;         // -1 : indicating that it has not been calculated yet
#     ifdef LOAD_BALANCE
      for (int v=0; v<4; v++)
      ParMassPos_Copy[v] = NULL;
#     else
      ParList_Copy = NULL;
#     endif

      for (int s=0; s<26; s++)
      {
         NPar_Escp   [s] = -1;   // -1 : indicating that it has not been calculated yet
         ParList_Escp[s] = NULL;
      }
#     endif

   } // METHOD : Activate



   //===================================================================================
   // Destructor  :  ~patch_t
   // Description :  Destructor of the structure "patch_t"
   //
   // Note        :  Deallocate flux and data arrays
   //===================================================================================
   ~patch_t()
   {

      fdelete();
      hdelete();
#     ifdef MHD
      edelete();
      mdelete();
#     endif
#     ifdef GRAVITY
      gdelete();
#     endif
#     ifdef DUAL_ENERGY
      sdelete();
#     endif

#     ifdef PARTICLE
      ddelete();

//    check: all particle-related variables/arrays should be deallocated already (except rho_ext when OPT__REUSE_MEMORY is on)
#     ifdef DEBUG_PARTICLE
      if ( ParList != NULL )              Aux_Error( ERROR_INFO, "ParList != NULL !!\n" );
#     ifdef LOAD_BALANCE
      for (int v=0; v<4; v++)
      if ( ParMassPos_Copy[v] != NULL )   Aux_Error( ERROR_INFO, "ParMassPos_Copy[%d] != NULL !!\n", v );
#     else
      if ( ParList_Copy != NULL )         Aux_Error( ERROR_INFO, "ParList_Copy != NULL !!\n" );
#     endif
      for (int s=0; s<26; s++)
      if ( ParList_Escp[s] != NULL )      Aux_Error( ERROR_INFO, "ParList_Escp[%d] != NULL !!\n", s );

      if ( NPar != 0 )                    Aux_Error( ERROR_INFO, "NPar = %d != 0 !!\n", NPar );
      if ( NParType[0] != 0 )             Aux_Error( ERROR_INFO, "NParType[0] = %d != 0 !!\n", NParType[0] );
      if ( NParType[1] != 0 )             Aux_Error( ERROR_INFO, "NParType[1] = %d != 0 !!\n", NParType[1] );
      if ( NPar_Copy != -1 )              Aux_Error( ERROR_INFO, "NPar_Copy = %d != -1 !!\n", NPar_Copy );
      for (int s=0; s<26; s++)
      if ( NPar_Escp[s] != -1 )           Aux_Error( ERROR_INFO, "NPar_Escp[%d] = %d != -1 !!\n", s, NPar_Escp[s] );
#     endif // #ifdef DEBUG_PARTICLE
#     endif // #ifdef PARTICLE

   } // METHOD : ~patch_t



   //===================================================================================
   // Method      :  fnew
   // Description :  Allocate flux[] in the given direction
   //
   // Note        :  flux[] will be initialized as zero
   //
   // Parameter   :  SibID    : Targeted sibling direction (0,1,2,3,4,5) <--> (-x,+x,-y,+y,-z,+z)
   //                AllocTmp : Allocate the temporary flux array flux_tmp[]
   //===================================================================================
   void fnew( const int SibID, const bool AllocTmp )
   {

#     ifdef GAMER_DEBUG
      if ( SibID < 0  ||  SibID > 5 )
         Aux_Error( ERROR_INFO, "incorrect SibID (%d) in fnew() !!\n", SibID );

      if ( flux[SibID] != NULL )
         Aux_Error( ERROR_INFO, "flux[%d] already exists !!\n", SibID );

      if ( AllocTmp  &&  flux_tmp[SibID] != NULL )
         Aux_Error( ERROR_INFO, "flux_tmp[%d] already exists !!\n", SibID );

#     ifdef BIT_REP_FLUX
      if ( flux_bitrep[SibID] != NULL )
         Aux_Error( ERROR_INFO, "flux_bitrep[%d] already exists !!\n", SibID );
#     endif
#     endif

      flux      [SibID]  = new real [NFLUX_TOTAL][PS1][PS1];
      if ( AllocTmp )
      flux_tmp  [SibID]  = new real [NFLUX_TOTAL][PS1][PS1];
#     ifdef BIT_REP_FLUX
      flux_bitrep[SibID] = new real [NFLUX_TOTAL][PS1][PS1];
#     endif

      for(int v=0; v<NFLUX_TOTAL; v++)
      for(int m=0; m<PS1; m++)
      for(int n=0; n<PS1; n++)
      {
         flux       [SibID][v][m][n] = 0.0;
         /*
//       no need to initialize flux_tmp[]
         if ( AllocTmp )
         flux_tmp   [SibID][v][m][n] = 0.0;
         */
#        ifdef BIT_REP_FLUX
         flux_bitrep[SibID][v][m][n] = 0.0;
#        endif
      }

   } // METHOD : fnew



   //===================================================================================
   // Method      :  fdelete
   // Description :  Deallocate flux[] along all directions
   //===================================================================================
   void fdelete()
   {

      for (int s=0; s<6; s++)
      {
         delete [] flux[s];
         flux[s] = NULL;

         delete [] flux_tmp[s];
         flux_tmp[s] = NULL;

#        ifdef BIT_REP_FLUX
         delete [] flux_bitrep[s];
         flux_bitrep[s] = NULL;
#        endif
      }

   } // METHOD : fdelete



#  ifdef MHD
   //===================================================================================
   // Method      :  enew
   // Description :  Allocate electric[] in the given direction
   //
   // Note        :  electric[] will be initialized as zero
   //
   // Parameter   :  SibID    : Target sibling direction (0-17)
   //                AllocTmp : Allocate the temporary electric array electric_tmp[]
   //===================================================================================
   void enew( const int SibID, const bool AllocTmp )
   {

#     ifdef GAMER_DEBUG
      if ( SibID < 0  ||  SibID > 17 )
         Aux_Error( ERROR_INFO, "incorrect SibID (%d) in fnew() !!\n", SibID );

      if ( electric[SibID] != NULL )
         Aux_Error( ERROR_INFO, "electric[%d] already exists !!\n", SibID );

      if ( AllocTmp  &&  electric_tmp[SibID] != NULL )
         Aux_Error( ERROR_INFO, "electric_tmp[%d] already exists !!\n", SibID );

#     ifdef BIT_REP_ELECTRIC
      if ( electric_bitrep[SibID] != NULL )
         Aux_Error( ERROR_INFO, "electric_bitrep[%d] already exists !!\n", SibID );
#     endif
#     endif

      const int Size = ( SibID < 6 ) ? NCOMP_ELE*PS1M1*PS1 : PS1;

      electric      [SibID]  = new real [Size];
      if ( AllocTmp )
      electric_tmp  [SibID]  = new real [Size];
#     ifdef BIT_REP_ELECTRIC
      electric_bitrep[SibID] = new real [Size];
#     endif

      for(int t=0; t<Size; t++)
      {
         electric       [SibID][t] = 0.0;
         /*
//       no need to initialize electric_tmp[]
         if ( AllocTmp )
         electric_tmp   [SibID][t] = 0.0;
         */
#        ifdef BIT_REP_ELECTRIC
         electric_bitrep[SibID][t] = 0.0;
#        endif
      }

   } // METHOD : enew



   //===================================================================================
   // Method      :  edelete
   // Description :  Deallocate electric[] along all directions
   //===================================================================================
   void edelete()
   {

      for (int s=0; s<18; s++)
      {
         delete [] electric[s];
         electric[s] = NULL;

         delete [] electric_tmp[s];
         electric_tmp[s] = NULL;

#        ifdef BIT_REP_ELECTRIC
         delete [] electric_bitrep[s];
         electric_bitrep[s] = NULL;
#        endif
      }

   } // METHOD : edelete
#  endif // #ifdef MHD



   //===================================================================================
   // Method      :  hnew
   // Description :  Allocate fluid[]
   //
   // Note        :  Do nothing if fluid[] has been allocated
   //===================================================================================
   void hnew()
   {

      if ( fluid == NULL )
      {
         fluid = new real [NCOMP_TOTAL][PS1][PS1][PS1];
         fluid[0][0][0][0] = (real)-1.0;  // arbitrarily initialized
      }

   } // METHOD : hnew



   //===================================================================================
   // Method      :  hdelete
   // Description :  Deallocate fluid[]
   //===================================================================================
   void hdelete()
   {

      delete [] fluid;
      fluid = NULL;

#     ifdef PARTICLE
      delete [] rho_ext;
      rho_ext = NULL;
#     endif

   } // METHOD : hdelete



#  ifdef MHD
   //===================================================================================
   // Method      :  mnew
   // Description :  Allocate magnetic[]
   //
   // Note        :  Do nothing if magnetic[] has been allocated
   //===================================================================================
   void mnew()
   {

      if ( magnetic == NULL )
      {
         magnetic = new real [NCOMP_MAG][ PS1P1*SQR(PS1) ];
         magnetic[0][0] = (real)-1.0;  // arbitrarily initialized
      }

   } // METHOD : mnew



   //===================================================================================
   // Method      :  mdelete
   // Description :  Deallocate magnetic[]
   //===================================================================================
   void mdelete()
   {

      delete [] magnetic;
      magnetic = NULL;

   } // METHOD : mdelete
#  endif // #ifdef MHD



#  ifdef GRAVITY
   //===================================================================================
   // Method      :  gnew
   // Description :  Allocate pot[] (and pot_ext[] for STORE_POT_GHOST)
   //
   // Note        :  Do nothing if pot[] (and pot_ext[] for STORE_POT_GHOST) has been allocated
   //===================================================================================
   void gnew()
   {

      if ( pot == NULL )      pot     = new real [PS1][PS1][PS1];

#     ifdef STORE_POT_GHOST
      if ( pot_ext == NULL )  pot_ext = new real [GRA_NXT][GRA_NXT][GRA_NXT];

//    always initialize pot_ext[] (even if pot_ext != NULL when calling this function) to indicate that this array
//    has NOT been properly set --> used by Poi_StorePotWithGhostZone()
      pot_ext[0][0][0] = POT_EXT_NEED_INIT;
#     endif

   } // METHOD : gnew



   //===================================================================================
   // Method      :  gdelete
   // Description :  Deallocate pot[] (and pot_ext[] for STORE_POT_GHOST)
   //===================================================================================
   void gdelete()
   {

      delete [] pot;
      pot = NULL;

#     ifdef STORE_POT_GHOST
      delete [] pot_ext;
      pot_ext = NULL;
#     endif

   } // METHOD : gdelete
#  endif // #ifdef GRAVITY



#  ifdef DUAL_ENERGY
   //===================================================================================
   // Method      :  snew
   // Description :  Allocate the dual-energy status array de_status[]
   //
   // Note        :  Do nothing if de_status[] has been allocated
   //===================================================================================
   void snew()
   {

      if ( de_status == NULL )
      {
         de_status = new char [PS1][PS1][PS1];
      }

   } // METHOD : snew



   //===================================================================================
   // Method      :  sdelete
   // Description :  Deallocate the dual-energy status array de_status[]
   //===================================================================================
   void sdelete()
   {

      delete [] de_status;
      de_status = NULL;

   } // METHOD : sdelete
#  endif // #ifdef DUAL_ENERGY



#  ifdef PARTICLE
   //===================================================================================
   // Method      :  dnew
   // Description :  Allocate rho_ext[]
   //
   // Note        :  Do nothing if rho_ext[] has been allocated
   //===================================================================================
   void dnew()
   {

      if ( rho_ext == NULL )  rho_ext = new real [RHOEXT_NXT][RHOEXT_NXT][RHOEXT_NXT];

//    always initialize rho_ext (even if rho_ext != NULL when calling this this function) to indicate that this array
//    has NOT been properly set --> used by Prepare_PatchData()
      rho_ext[0][0][0] = RHO_EXT_NEED_INIT;

   } // METHOD : dnew



   //===================================================================================
   // Method      :  ddelete
   // Description :  Deallocate rho_ext[]
   //===================================================================================
   void ddelete()
   {

      delete [] rho_ext;
      rho_ext = NULL;

   } // METHOD : ddelete



   //===================================================================================
   // Method      :  AddParticle
   // Description :  Add new particles into the particle list
   //
   // Note        :  1. Particles in the list must already lie inside the target patch
   //                   --> It will be checked in the debug mode
   //                2. If ParList already exists, new particles will be attached to the existing particle list
   //                3. Also update the total number of particles at the target level
   //
   // Parameter   :  NNew    : Number of new particles to be added
   //                NewList : List storing the indices of new particles
   //                NPar_Lv : Pointer to amr->Par->NPar_Lv[TargetLv]
   //                ParPos  : Particle position list             (debug only)
   //                NParTot : Total number of existing particles (debug only)
   //                Comment : Message for the debug mode         (debug only)
   //===================================================================================
#  ifdef DEBUG_PARTICLE
   void AddParticle( const int NNew, const long *NewList, long *NPar_Lv,
                     const real **ParPos, const long NParTot, const char *Comment )
#  else
   void AddParticle( const int NNew, const long *NewList, long *NPar_Lv )
#  endif
   {

//    check
#     ifdef DEBUG_PARTICLE
//    check 1: no strange settings
      if ( NewList == NULL  &&  NNew != 0 )  Aux_Error( ERROR_INFO, "\"%s\": NewList == NULL !!\n", Comment );
      if ( NNew < 0 )                        Aux_Error( ERROR_INFO, "\"%s\": NNew (%d) < 0 !!\n",   Comment, NNew );
      if ( NPar < 0 )                        Aux_Error( ERROR_INFO, "\"%s\": NPar (%d) < 0 !!\n",   Comment, NPar );
      if ( NParType[0] < 0 )                 Aux_Error( ERROR_INFO, "\"%s\": NParType[0] (%d) < 0 !!\n", Comment, NParType[0] );
      if ( NParType[1] < 0 )                 Aux_Error( ERROR_INFO, "\"%s\": NParType[1] (%d) < 0 !!\n", Comment, NParType[1] );

//    check 2: particle indices in NewList are NEW
      for (int q=0; q<NPar; q++)
      for (int p=0; p<NNew; p++)
      {
         if ( ParList[q] == NewList[p] )
            Aux_Error( ERROR_INFO, "\"%s\": repeated particle index (NewList[%d] = %ld) !!\n", Comment, p, NewList[p] );
      }

//    check 3: no duplicate particle index in NewList
      for (int q=0; q<NNew; q++)
      for (int p=0; p<NNew; p++)
      {
         if ( p != q  &&  NewList[q] == NewList[p] )
            Aux_Error( ERROR_INFO, "\"%s\": duplicate particle index (NewList[%d] = NewList[%d] = %ld) !!\n",
                       Comment, p, q, NewList[p] );
      }

//    check 4: whether the new particles indeed lie inside the target patch
      for (int p=0; p<NNew; p++)
      {
         const long ParID = NewList[p];

         for (int d=0; d<3; d++)
         {
            if ( ParPos[d][ParID] != ParPos[d][ParID]  ||  ParPos[d][ParID] < EdgeL[d]  ||  ParPos[d][ParID] >= EdgeR[d] )
               Aux_Error( ERROR_INFO, "\"%s\": wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d][%d] = %13.6e) !!\n",
                          Comment, EdgeL[d], EdgeR[d], d, ParID, ParPos[d][ParID] );
         }
      }

//    check 5: new particle IDs do not exceed the total number of particles
      for (int p=0; p<NNew; p++)
      {
         const long ParID = NewList[p];

         if ( ParID >= NParTot )
            Aux_Error( ERROR_INFO, "\"%s\": Particle ID (%ld) >= Total number of active + inactive particles (%ld) !!\n",
                       Comment, ParID, NParTot );
      }

      if ( NPar_Lv == NULL)   Aux_Error( ERROR_INFO, "NPar_Lv == NULL !!\n" );
#     endif // #ifdef DEBUG_PARTICLE


//    nothing to do if NNew == 0
      if ( NNew == 0 )  return;


//    allocate enough memory for the particle list array
      const int NPar_New = NPar + NNew;
      if ( NPar_New > ParListSize )
      {
         ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*NPar_New );
         ParList     = (long*)realloc( ParList, ParListSize*sizeof(long) );
      }

//    record the new particle indices
      for (int p=0; p<NNew; p++)    ParList[ NPar + p ] = NewList[p];

//    update the particle number
      NPar = NPar_New;

//    update the particle type number

      for (int p=0; p<NPar; p++) NParType[(int)ParList[p].Type]++

//    update the particle number at the target level
      *NPar_Lv += NNew;

   } // METHOD : AddParticle



   //===================================================================================
   // Method      :  RemoveParticle
   // Description :  Remove particles from the particle list
   //
   // Note        :  1. The numbers stored in RemoveList are "array indices of ParList"
   //                   instead of "particle indices stored in ParList"
   //                2. The numbers stored in RemoveList must be in ascending numerical order
   //                3. Also update the total number of particles at the target level
   //
   // Parameter   :  NRemove    : Number of particles to be removed
   //                RemoveList : List storing the array indices of "ParList" to be removed
   //                             **array indices, NOT particle indices**
   //                NPar_Lv    : Pointer to amr->Par->NPar_Lv[TargetLv]
   //                             --> NPar_Lv will NOT be updated if it's NULL
   //                             --> This is useful for OpenMP, for which different threads
   //                                 may tend to update NPar_Lv at the same time
   //                             --> So for NPar_Lv == NULL, one must update NPar_Lv later manually
   //                RemoveAll  : true --> remove all particle in this patch
   //===================================================================================
   void RemoveParticle( const int NRemove, const int *RemoveList, long *NPar_Lv, const bool RemoveAll )
   {

//    removing all particles is easy
      if ( RemoveAll )
      {
//       update the particle number at the target level
         if ( NPar_Lv != NULL )
         {
            *NPar_Lv -= NPar;

#           ifdef DEBUG_PARTICLE
            if ( *NPar_Lv < 0 )  Aux_Error( ERROR_INFO, "NPar_Lv = %ld < 0 !!\n", *NPar_Lv );
#           endif
         }


//       remove all particles
         NPar        = 0;
         ParListSize = 0;
         NParType[0] = 0;
         NParType[1] = 0;

         if ( ParList != NULL )
         {
            free( ParList );
            ParList = NULL;
         }

         return;
      } // if ( RemoveAll )


//    check
#     ifdef DEBUG_PARTICLE
//    check 1: no strange settings
      if ( RemoveList == NULL  &&  NRemove != 0 )  Aux_Error( ERROR_INFO, "RemoveList == NULL !!\n" );
      if ( NRemove < 0 )                           Aux_Error( ERROR_INFO, "NRemove (%d) < 0 !!\n", NRemove );
      if ( NPar < 0 )                              Aux_Error( ERROR_INFO, "NPar (%d) < 0 !!\n", NPar );
      if ( NParType[0] < 0 )                       Aux_Error( ERROR_INFO, "NParType[0] (%d) < 0 !!\n", NParType[0] );
      if ( NParType[1] < 0 )                       Aux_Error( ERROR_INFO, "NParType[1] (%d) < 0 !!\n", NParType[1] );
      if ( NRemove > NPar )                        Aux_Error( ERROR_INFO, "NRemove (%d) > NPar (%d) !!\n", NRemove, NPar );

//    check 2: indices in RemoveList do not exceed the number of existing particles
      for (int p=0; p<NRemove; p++)
      {
         if ( RemoveList[p] >= NPar  ||  RemoveList[p] < 0 )
            Aux_Error( ERROR_INFO, "incorrect RemoveList[%d]=%d, (NPar=%d) !!\n", p, RemoveList[p], NPar );
      }

//    check 3: no duplicate index in RemoveList
      for (int q=0; q<NRemove; q++)
      for (int p=0; p<NRemove; p++)
      {
         if ( p != q  &&  RemoveList[q] == RemoveList[p] )
            Aux_Error( ERROR_INFO, "duplicate index (RemoveList[%d] = RemoveList[%d] = %d) !!\n",
                       p, q, RemoveList[p] );
      }

//    check 4: RemoveList is in ascending numerical order
      for (int p=1; p<NRemove; p++)
      {
         if ( RemoveList[p] <= RemoveList[p-1] )
            Aux_Error( ERROR_INFO, "RemoveList is NOT in ascending numerical order ([%d]=%d <= [%d]=%d) !!\n",
                       p, RemoveList[p], p-1, RemoveList[p-1] );
      }
#     endif // #ifdef DEBUG_PARTICLE


//    nothing to do if NRemove == 0
      if ( NRemove == 0 )  return;


//    remove particles
      int LastP;
      for (int p=NRemove-1; p>=0; p--)
      {
         LastP = NPar - 1;

         NParType[(int)ParList[RemoveList[p]].Type] --;

         if ( RemoveList[p] != LastP )   ParList[ RemoveList[p] ] = ParList[ LastP ];

         NPar --;
      }


//    release some memory if we waste too much
      const int ParListSize_Min = (int)floor( PARLIST_REDUCE_FACTOR*ParListSize );

      if ( NPar <= ParListSize_Min )
      {
         ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*NPar );
         ParList     = (long*)realloc( ParList, ParListSize*sizeof(long) );
      }


//    update the particle number at the target level
      if ( NPar_Lv != NULL )
      {
         *NPar_Lv -= NRemove;

#        ifdef DEBUG_PARTICLE
         if ( *NPar_Lv < 0 )  Aux_Error( ERROR_INFO, "NPar_Lv = %ld < 0 !!\n", *NPar_Lv );
#        endif
      }

   } // METHOD : RemoveParticle
#  endif // #ifdef PARTICLE


}; // struct patch_t



#endif // #ifndef __PATCH_H__
