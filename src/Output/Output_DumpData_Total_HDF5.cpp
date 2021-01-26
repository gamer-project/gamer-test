#ifdef SUPPORT_HDF5

#include "GAMER.h"
#include "HDF5_Typedef.h"
#include <ctime>

void FillIn_KeyInfo  (   KeyInfo_t &KeyInfo   );
void FillIn_Makefile (  Makefile_t &Makefile  );
void FillIn_SymConst (  SymConst_t &SymConst  );
void FillIn_InputPara( InputPara_t &InputPara );

static void GetCompound_KeyInfo  ( hid_t &H5_TypeID );
static void GetCompound_Makefile ( hid_t &H5_TypeID );
static void GetCompound_SymConst ( hid_t &H5_TypeID );
static void GetCompound_InputPara( hid_t &H5_TypeID );



/*======================================================================================================
Data structure:
/ -> |
     | -> Info group     -> | -> InputPara dset (compound)
     |                      | -> KeyInfo   dset (compound)
     |                      | -> Makefile  dset (compound)
     |                      | -> SymConst  dset (compound)
     |
     | -> Tree group     -> | -> Corner  dset -> Cvt2Phy attrs
     |                      | -> LBIdx   dset
     |                      | -> Father  dset
     |                      | -> Son     dset
     |                      | -> Sibling dset
     |                      | -> NPar    dset
     |
     | -> GridData group -> | -> Dens dset
     |                      | -> ...
     |                      | -> ...
     |
     | -> Particle group -> | -> ParMass dset
                            | -> ...
                            | -> ...
======================================================================================================*/



/*======================================================================================================
h5py usage (with python 2):
1. Load file: f=h5py.File("Data_000000", "r")
2. Shows the names of all groups: list(f) or f.keys()
3. Show the names of all attributes: list(f['Tree']['Corner'].attrs) or f['Tree']['Corner'].attrs.keys()
4. Show a specific attribute: f['Tree']['Corner'].attrs['Cvt2Phy']
5. Show all variables in a compound variable: f['Info']['KeyInfo'].dtype
6. Show the value of a dataset: f['Tree']['Corner'].value
7. Show density of a patch with a global ID (GID) 1234: f['Data']['Dens'][1234]
8. Show density at a specific cell [1][2][3]: f['Data']['Dens'][1234][1][2][3]
======================================================================================================*/



/*======================================================================================================
Procedure for outputting new variables:
1. Add the new variable into one of the data structures (XXX) defined in "HDF5_Typedef.h"
2. Edit "GetCompound_XXX" to insert the new variables into the compound datatype
3. Edit "FillIn_XXX" to fill in the new variables
4. Edit "Check_XXX" in "Init_ByRestart_HDF5.cpp" to load and compare the new variables
5. Update FormatVersion
======================================================================================================*/




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Total_HDF5 (FormatVersion = 2429)
// Description :  Output all simulation data in the HDF5 format, which can be used as a restart file
//                or loaded by YT
//
// Note        :  1. Please refer to the "Data structure" described on the top of this file
//                2. Patch IDs stored in the HDF5 output are always GID (global identification) instead of
//                   PID (patch identification). Unlike PID, which always starts from 0 at different ranks
//                   and different levels, GID is unique among all patches at all ranks and all levels
//                   --> Each patch has an unique GID
//                3. Both "Father, Son, and Sibling[26]" are GID instead of PID
//                4. Currently we always use HDF5 NATVIE datatypes for both memory and dataset
//                5. All arrays in the "Tree" group (e.g., Corner, LBIdx, ...) and the "GridData" group (e.g., Dens, MomX, ...)
//                   have been sorted according to GID
//                   --> Moreover, currently we store all patches at the same level together
//                   --> A higher-level patch always has a larger GID than a lower-level patch
//                5. LBIdx dataset stores the LB_Idx of all patches sorted by their GIDs
//                   --> This list will be created even when LOAD_BALANCE is not turned on
//                       so that a serial output can be loaded by a parallel run easily
//                6. All C structures (e.g., "KeyInfo, SymConst, ..." are stored as a single (scalar)
//                   compound datetype
//                   --> Better update h5py to version >= 2.3.0 to properly read it in python
//                   --> Corresponding C structures are defined in "HDF5_Typedef.h"
//                7. It seems that h5py still have problem for "modifying" the loaded data. But reading data is fine.
//                8. The "H5T_GAMER_REAL" datatype will be mapped to "H5T_NATIVE_DOUBLE / H5T_NATIVE_FLOAT" if
//                   FLOAT8 is on / off
//                9. It is found that in the parallel environment each rank must try to "synchronize" the HDF5 file
//                   before opening the existed file and add data
//                   --> To achieve that, always invoke "SyncHDF5File" before calling "H5Fopen"
//                   --> "SyncHDF5File" is defined in "HDF5_Typedef.h", which simply opens the file
//                       with the appending mode and then closes it immediately
//                10. With PARTICLE on, two additional particle information will be recorded:
//                    --> "NPar" dataset under "Tree" records the number of active particles in all patches
//                         sorted by their GIDs
//                    --> "Particle" dataset under "/" stores all particle attributes (e.g., mass, velocity, ...)
//                        --> Currently we store different attributes in separate datasets
//                        --> Particles are stored in the order of their associated GIDs as well, but the order of
//                            particles in the same patch is not specified
//
// Parameter   :  FileName : Name of the output file
//
// Revision    :  2210 : 2016/10/03 --> output HUBBLE0, OPT__UNIT, UNIT_L/M/T/V/D/E, MOLECULAR_WEIGHT
//                2216 : 2016/11/27 --> output OPT__FLAG_LOHNER_TEMP
//                2217 : 2017/01/25 --> output RESTART_LOAD_NRANK, set CodeVersion to "gamer"
//                2218 : 2017/01/28 --> output OPT__FLAG_VORTICITY and the corresponding flag table
//                2220 : 2017/02/14 --> output passive grid and particle variables
//                2221 : 2017/02/20 --> output TINY_NUMBER and HUGE_NUMBER
//                2222 : 2017/02/20 --> output OPT__NORMALIZE_PASSIVE
//                2223 : 2017/02/22 --> output NormalizePassive_NVar and NormalizePassive_VarIdx
//                2224 : 2017/02/25 --> output OPT__CK_NORMALIZE_PASSIVE
//                2225 : 2017/03/01 --> output LB_Par_Weight, rename LB_Input__WLI_Max as LB_WLI_Max
//                2226 : 2017/03/03 --> output Opt__RecordLoadBalance
//                2227 : 2017/03/21 --> output PassiveFieldName_Grid and PassiveFieldName_Par
//                2228 : 2017/03/21 --> output NCOMP_FLUID, NCOMP_PASSIVE, and PAR_NPASSIVE in KeyInfo_t
//                2229 : 2017/04/06 --> output DUAL_ENERGY and DUAL_ENERGY_SWITCH
//                2230 : 2017/05/08 --> output OPT__FLAG_PAR_MASS_CELL and FlagTable_ParMassCell
//                2231 : 2017/05/08 --> output OPT__FLAG_JEANS and FlagTable_Jeans
//                2232 : 2017/06/13 --> output TESTPROB_ID
//                2233 : 2017/06/13 --> rename Opt__Output_TestError as Opt__Output_User
//                2234 : 2017/06/25 --> output Grackle variables
//                2235 : 2017/07/11 --> replace IndividualDt and Opt__AdaptiveDt by Opt__DtLevel
//                2236 : 2017/07/16 --> output DT_FLU_BLOCK_SIZE and DT_FLU_USE_SHUFFLE
//                2237 : 2017/07/17 --> output DT__FLEXIBLE_RANGE
//                2238 : 2017/07/19 --> output DT_GRA_BLOCK_SIZE_Z and DT_FLU_USE_SHUFFLE
//                2239 : 2017/07/25 --> output JEANS_MIN_PRES, JEANS_MIN_PRES_LEVEL, JEANS_MIN_PRES_NCELL
//                2240 : 2017/07/26 --> output AUTO_REDUCE_DT*
//                2241 : 2017/08/05 --> output FLAG_BUFFER_SIZE_MAXM1_LV
//                2242 : 2017/08/09 --> output DT__SYNC_PARENT_LV and DT__SYNC_CHILDREN_LV; remove DT__FLEXIBLE_RANGE
//                2250 : 2017/08/09 --> output dTime_AllLv[]
//                2251 : 2017/08/12 --> output OPT__RESTART_RESET
//                2252 : 2017/08/17 --> output GRACKLE_PE_HEATING and GRACKLE_PE_HEATING_RATE
//                2253 : 2017/08/27 --> output STAR_FORMATION
//                2254 : 2017/08/29 --> output star formation parameters
//                2255 : 2017/09/01 --> output SF_CREATION_STAR_DET_RANDOM
//                2256 : 2017/09/10 --> output FLAG_BUFFER_SIZE_MAXM2_LV
//                2257 : 2017/09/17 --> output OPT__OPTIMIZE_AGGRESSIVE
//                2258 : 2017/09/21 --> output OPT__MINIMIZE_MPI_BARRIER
//                2259 : 2017/10/10 --> output OPT__INIT_GRID_WITH_OMP
//                2261 : 2017/12/05 --> no longer define INTEL
//                2262 : 2017/12/27 --> rename all UM variables
//                2263 : 2017/12/27 --> remove OPT__RESTART_HEADER
//                2264 : 2018/02/28 --> add RANDOM_NUMBER
//                2265 : 2018/04/02 --> add OPT__NO_FLAG_NEAR_BOUNDARY
//                2266 : 2018/05/11 --> add OPT__UM_IC_LOAD_NRANK
//                2300 : 2018/07/15 --> replace PAR_NVAR and PAR_NPASSIVE by PAR_NATT_STORED and PAR_NATT_USER;
//                                      use the new infrastructure for adding user-defined grid fields and
//                                      particle attributes
//                                      --> imcompatible with version 2266 for the data with user-defined grid fields
//                                          and particle attributes as their labels may have changed
//                2301 : 2018/07/24 --> add OPT__UM_IC_FORMAT, PAR_IC_FORMAT, and PAR_IC_MASS
//                2302 : 2018/07/24 --> replace GRACKLE_MODE by GRACKLE_ACTIVATE
//                2303 : 2018/10/04 --> set "CodeVersion" to VERSION defined in Macro.h
//                2304 : 2018/12/10 --> remove EP_Coeff that no longer exists
//                2305 : 2018/12/15 --> remove variables related to the WAF scheme
//                2306 : 2018/12/25 --> replace DT_GRA_BLOCK_SIZE_Z by DT_GRA_BLOCK_SIZE
//                2307 : 2018/12/27 --> replace GRA_BLOCK_SIZE_Z by GRA_BLOCK_SIZE
//                2308 : 2019/03/14 --> add OPT__RECORD_NOTE and OPT__RECORD_UNPHY
//                2309 : 2019/03/27 --> add OPT__FIXUP_ELECTRIC
//                2310 : 2019/04/20 --> add OPT__CK_INTERFACE_B
//                2311 : 2019/05/22 --> add OPT__CK_DIVERGENCE_B
//                2312 : 2019/05/31 --> add OPT__GRAVITY_EXTRA_MASS
//                2400 : 2019/06/08 --> output magnetic field for MHD
//                2401 : 2019/06/30 --> output OPT__FLAG_CURRENT and FlagTable_Current
//                2402 : 2019/07/17 --> replace USG_GhostSize by USG_GhostSizeF and USG_GhostSizeG
//                2403 : 2019/09/20 --> add BIT_REP_FLUX and BIT_REP_ELECTRIC defined in CUFLU.h
//                2404 : 2019/10/16 --> add DT__MAX
//                2405 : 2019/12/29 --> output GRACKLE_THREE_BODY_RATE, GRACKLE_CIE_COOLING, GRACKLE_H2_OPA_APPROX
//                2406 : 2020/02/26 --> output EOS
//                2407 : 2020/02/27 --> output MIN_EINT
//                2408 : 2020/05/10 --> output EOS_NAUX_MAX
//                2409 : 2020/07/05 --> output HLLC_WAVESPEED
//                2410 : 2020/07/31 --> output HLLE_WAVESPEED
//                2411 : 2020/08/11 --> output LR_EINT
//                2412 : 2020/08/17 --> output FLU_NIN_T
//                2413 : 2020/08/23 --> output HLLD_WAVESPEED
//                2414 : 2020/09/06 --> output INT_OPP_SIGN_0TH_ORDER
//                2415 : 2020/09/08 --> output OPT__LAST_RESORT_FLOOR
//                2416 : 2020/09/08 --> output BAROTROPIC_EOS
//                2417 : 2020/09/09 --> output ISO_TEMP
//                2418 : 2020/09/21 --> replace (OPT__GRAVITY_TYPE, OPT__EXTERNAL_POT) by
//                                      (OPT__SELF_GRAVITY, OPT__EXT_ACC, OPT__EXT_POT)
//                2419 : 2020/09/25 --> output EXTPOT_BLOCK_SIZE
//                2420 : 2020/10/12 --> output OPT__FLAG_USER_NUM and use variable-length datatype for FlagTable_User
//                2421 : 2020/10/26 --> output COSMIC_RAY
//                2422 : 2020/10/29 --> output the parameters of external potential table
//                2423 : 2020/11/01 --> output EOS_NTABLE_MAX
//                2424 : 2020/12/22 --> output SRC_USER
//                2425 : 2020/12/24 --> output SRC_DELEPTONIZATION
//                2426 : 2020/12/24 --> output FLU_NIN_S, FLU_NOUT_S, and SRC_GPU_NPGROUP
//                2427 : 2020/12/26 --> output SRC_BLOCK_SIZE and SRC_GHOST_SIZE
//                2428 : 2020/12/27 --> output SRC_NAUX_DLEP and SRC_NAUX_USER
//                2429 : 2021/01/26 --> output SRC_DLEP_PROF_NVAR and SRC_DLEP_PROF_NBINMAX
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Total_HDF5( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the target file already exists
   if ( Aux_CheckFileExist(FileName)  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// 1. gather the number of patches at different MPI ranks and set the corresponding GID offset
   int (*NPatchAllRank)[NLEVEL] = new int [MPI_NRank][NLEVEL];
   int NPatchLocal[NLEVEL], NPatchAllLv=0, GID_Offset[NLEVEL], GID_LvStart[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)  NPatchLocal[lv] = amr->NPatchComma[lv][1];

   MPI_Allgather( NPatchLocal, NLEVEL, MPI_INT, NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      GID_Offset[lv] += NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   GID_Offset[lv] += NPatchTotal[FaLv];

      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }



// 2. prepare all HDF5 variables
   hsize_t H5_SetDims_LBIdx, H5_SetDims_Cr[2], H5_SetDims_Fa, H5_SetDims_Son, H5_SetDims_Sib[2], H5_SetDims_Field[4];
   hsize_t H5_MemDims_Field[4], H5_Count_Field[4], H5_Offset_Field[4];
   hid_t   H5_MemID_Field;
   hid_t   H5_FileID, H5_GroupID_Info, H5_GroupID_Tree, H5_GroupID_GridData;
   hid_t   H5_SetID_LBIdx, H5_SetID_Cr, H5_SetID_Fa, H5_SetID_Son, H5_SetID_Sib, H5_SetID_Field;
   hid_t   H5_SetID_KeyInfo, H5_SetID_Makefile, H5_SetID_SymConst, H5_SetID_InputPara;
   hid_t   H5_SpaceID_Scalar, H5_SpaceID_LBIdx, H5_SpaceID_Cr, H5_SpaceID_Fa, H5_SpaceID_Son, H5_SpaceID_Sib, H5_SpaceID_Field;
   hid_t   H5_TypeID_Com_KeyInfo, H5_TypeID_Com_Makefile, H5_TypeID_Com_SymConst, H5_TypeID_Com_InputPara;
   hid_t   H5_DataCreatePropList;
   hid_t   H5_AttID_Cvt2Phy;
   herr_t  H5_Status;
#  ifdef PARTICLE
   hsize_t H5_SetDims_NPar, H5_SetDims_ParData[1], H5_MemDims_ParData[1],  H5_Count_ParData[1], H5_Offset_ParData[1];
   hid_t   H5_SetID_NPar, H5_SpaceID_NPar, H5_SpaceID_ParData, H5_GroupID_Particle, H5_SetID_ParData, H5_MemID_ParData;
#  endif
#  ifdef MHD
   hsize_t H5_SetDims_FCMag[4], H5_MemDims_FCMag[4], H5_Count_FCMag[4], H5_Offset_FCMag[4];
   hid_t   H5_MemID_FCMag, H5_SetID_FCMag, H5_SpaceID_FCMag[NCOMP_MAG];
#  endif

// 2-1. do NOT write fill values to any dataset for higher I/O performance
   H5_DataCreatePropList = H5Pcreate( H5P_DATASET_CREATE );
   H5_Status             = H5Pset_fill_time( H5_DataCreatePropList, H5D_FILL_TIME_NEVER );

// 2-2. create the "compound" datatype
   GetCompound_KeyInfo  ( H5_TypeID_Com_KeyInfo   );
   GetCompound_Makefile ( H5_TypeID_Com_Makefile  );
   GetCompound_SymConst ( H5_TypeID_Com_SymConst  );
   GetCompound_InputPara( H5_TypeID_Com_InputPara );

// 2-3. create the "scalar" dataspace
   H5_SpaceID_Scalar = H5Screate( H5S_SCALAR );



// 3. output the simulation information
   if ( MPI_Rank == 0 )
   {
//    3-1. collect all information to be recorded
      KeyInfo_t   KeyInfo;
      Makefile_t  Makefile;
      SymConst_t  SymConst;
      InputPara_t InputPara;

      FillIn_KeyInfo  ( KeyInfo   );
      FillIn_Makefile ( Makefile  );
      FillIn_SymConst ( SymConst  );
      FillIn_InputPara( InputPara );


//    3-2. create the HDF5 file (overwrite the existing file)
      H5_FileID = H5Fcreate( FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to create the HDF5 file \"%s\" !!\n", FileName );


//    3-3. write the simulation info (note: dataset doesn't support VL datatype when the fill value is not defined)
      H5_GroupID_Info = H5Gcreate( H5_FileID, "Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Info < 0 )    Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Info" );

//    3-3-1. KeyInfo
      H5_SetID_KeyInfo   = H5Dcreate( H5_GroupID_Info, "KeyInfo", H5_TypeID_Com_KeyInfo, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_KeyInfo < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "KeyInfo" );
      H5_Status          = H5Dwrite( H5_SetID_KeyInfo, H5_TypeID_Com_KeyInfo, H5S_ALL, H5S_ALL, H5P_DEFAULT, &KeyInfo );
      H5_Status          = H5Dclose( H5_SetID_KeyInfo );

//    3-3-2. Makefile
      H5_SetID_Makefile  = H5Dcreate( H5_GroupID_Info, "Makefile", H5_TypeID_Com_Makefile, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_Makefile < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Makefile" );
      H5_Status          = H5Dwrite( H5_SetID_Makefile, H5_TypeID_Com_Makefile, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Makefile );
      H5_Status          = H5Dclose( H5_SetID_Makefile );

//    3-3-3. SymConst
      H5_SetID_SymConst  = H5Dcreate( H5_GroupID_Info, "SymConst", H5_TypeID_Com_SymConst, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_SymConst < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "SymConst" );
      H5_Status          = H5Dwrite( H5_SetID_SymConst, H5_TypeID_Com_SymConst, H5S_ALL, H5S_ALL, H5P_DEFAULT, &SymConst );
      H5_Status          = H5Dclose( H5_SetID_SymConst );

//    3-3-4. InputPara
      H5_SetID_InputPara = H5Dcreate( H5_GroupID_Info, "InputPara", H5_TypeID_Com_InputPara, H5_SpaceID_Scalar,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_SetID_InputPara < 0 ) Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "InputPara" );
      H5_Status          = H5Dwrite( H5_SetID_InputPara, H5_TypeID_Com_InputPara, H5S_ALL, H5S_ALL, H5P_DEFAULT, &InputPara );
      H5_Status          = H5Dclose( H5_SetID_InputPara );

      H5_Status = H5Gclose( H5_GroupID_Info );
      H5_Status = H5Fclose( H5_FileID );

//    3-4. free memory
      for (int lv=0; lv<NLEVEL-1; lv++)   free( InputPara.FlagTable_User[lv].p );
   } // if ( MPI_Rank == 0 )



// 4. output the AMR tree structure (father, son, sibling, LBIdx, corner, and the number of particles --> sorted by GID)
   long *LBIdxList_Local[NLEVEL], *LBIdxList_AllLv;
   int  (*CrList_Local[NLEVEL])[3], (*CrList_AllLv)[3];
   int  *FaList_Local[NLEVEL], *FaList_AllLv;
   int  *SonList_Local[NLEVEL], *SonList_AllLv;
   int  (*SibList_Local[NLEVEL])[26], (*SibList_AllLv)[26];
#  ifdef PARTICLE
   int  *NParList_Local[NLEVEL], *NParList_AllLv;
#  endif

   long *LBIdxList_Sort[NLEVEL];
   int  *LBIdxList_Sort_IdxTable[NLEVEL];

   int   MyGID, FaPID, FaGID, FaLv, SonPID, SonGID, SonLv, SibPID, SibGID, MatchIdx;
   long  FaLBIdx, SonLBIdx, SibLBIdx;
   int  *SonCr=NULL, *SibCr=NULL;

   int   RecvCount_LBIdx[MPI_NRank], RecvDisp_LBIdx[MPI_NRank], RecvCount_Cr[MPI_NRank], RecvDisp_Cr[MPI_NRank];
   int   RecvCount_Fa[MPI_NRank], RecvDisp_Fa[MPI_NRank], RecvCount_Son[MPI_NRank], RecvDisp_Son[MPI_NRank];
   int   RecvCount_Sib[MPI_NRank], RecvDisp_Sib[MPI_NRank];
#  ifdef PARTICLE
   int   RecvCount_NPar[MPI_NRank], RecvDisp_NPar[MPI_NRank];
#  endif

// 4-1. allocate lists
   if ( MPI_Rank == 0 )
   {
      LBIdxList_AllLv = new long [ NPatchAllLv ];
      CrList_AllLv    = new int  [ NPatchAllLv ][3];
      FaList_AllLv    = new int  [ NPatchAllLv ];
      SonList_AllLv   = new int  [ NPatchAllLv ];
      SibList_AllLv   = new int  [ NPatchAllLv ][26];
#     ifdef PARTICLE
      NParList_AllLv  = new int  [ NPatchAllLv ];
#     endif
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long [ amr->NPatchComma[lv][1] ];
      CrList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ][3];
      FaList_Local           [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SonList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ];
      SibList_Local          [lv] = new int  [ amr->NPatchComma[lv][1] ][26];
#     ifdef PARTICLE
      NParList_Local         [lv] = new int  [ amr->NPatchComma[lv][1] ];
#     endif

      LBIdxList_Sort         [lv] = new long [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int  [ NPatchTotal[lv] ];
   }


// 4-2. collect and sort LBIdx from all ranks
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_LBIdx[r] = NPatchAllRank[r][lv];
         RecvDisp_LBIdx [r] = ( r == 0 ) ? 0 : RecvDisp_LBIdx[r-1] + RecvCount_LBIdx[r-1];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

//    all ranks need to get LBIdxList_Sort since we will use it to calculate GID
      MPI_Allgatherv( LBIdxList_Local[lv], amr->NPatchComma[lv][1], MPI_LONG,
                      LBIdxList_Sort[lv], RecvCount_LBIdx, RecvDisp_LBIdx, MPI_LONG,
                      MPI_COMM_WORLD );
   } // for (int lv=0; lv<NLEVEL; lv++)

// store in the AllLv array BEFORE sorting
   if ( MPI_Rank == 0 )
   {
      MyGID = 0;

      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<NPatchTotal[lv]; PID++)
         LBIdxList_AllLv[ MyGID++ ] = LBIdxList_Sort[lv][PID];
   }

// sort list and get the corresponding index table (for calculating GID later)
   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_Sort[lv], LBIdxList_Sort_IdxTable[lv] );


// 4-3. store the local tree
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       4-3-1. LBIdx (set already)
//       LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;


//       4-3-2. corner
         for (int d=0; d<3; d++)
         CrList_Local[lv][PID][d] = amr->patch[0][lv][PID]->corner[d];


//       4-3-3. father GID
         FaPID = amr->patch[0][lv][PID]->father;
         FaLv  = lv - 1;

//       no father (only possible for the root patches)
         if ( FaPID < 0 )
         {
#           ifdef DEBUG_HDF5
            if ( lv != 0 )       Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 !!\n", lv, PID, FaPID );
            if ( FaPID != -1 )   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 but != -1 !!\n", lv, PID, FaPID );
#           endif

            FaGID = FaPID;
         }

//       father patch is a real patch
         else if ( FaPID < amr->NPatchComma[FaLv][1] )
            FaGID = FaPID + GID_Offset[FaLv];

//       father patch is a buffer patch (only possible in LOAD_BALANCE)
         else // (FaPID >= amr->NPatchComma[FaLv][1] )
         {
#           ifdef DEBUG_HDF5
#           ifndef LOAD_BALANCE
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= NRealFaPatch %d (only possible in LOAD_BALANCE) !!\n",
                       lv, PID, FaPID, amr->NPatchComma[FaLv][1] );
#           endif

            if ( FaPID >= amr->num[FaLv] )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= total number of patches %d !!\n",
                       lv, PID, FaPID, amr->num[FaLv] );
#           endif // DEBUG_HDF5

            FaLBIdx = amr->patch[0][FaLv][FaPID]->LB_Idx;

            Mis_Matching_int( NPatchTotal[FaLv], LBIdxList_Sort[FaLv], 1, &FaLBIdx, &MatchIdx );

#           ifdef DEBUG_HDF5
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d, FaLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, FaPID, FaLBIdx );
#           endif

            FaGID = LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + GID_LvStart[FaLv];
         } // if ( FaPID >= amr->NPatchComma[FaLv][1] )

         FaList_Local[lv][PID] = FaGID;


//       4-3-4. son GID
         SonPID = amr->patch[0][lv][PID]->son;
         SonLv  = lv + 1;

//       no son (must check this first since SonLv may be out of range --> == NLEVEL)
         if      ( SonPID == -1 )
            SonGID = SonPID;

//       son patch is a real patch at home
         else if ( SonPID >= 0  &&  SonPID < amr->NPatchComma[SonLv][1] )
            SonGID = SonPID + GID_Offset[SonLv];

//       son patch lives abroad (only possible in LOAD_BALANCE)
         else if ( SonPID < -1 )
         {
#           ifdef DEBUG_HDF5
#           ifdef LOAD_BALANCE
            const int SonRank = SON_OFFSET_LB - SonPID;
            if ( SonRank < 0  ||  SonRank == MPI_Rank  ||  SonRank >= MPI_NRank )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, incorrect SonRank %d (MyRank %d, NRank %d) !!\n",
                       lv, PID, SonPID, SonRank, MPI_Rank, MPI_NRank );
#           else
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d < -1 (only possible in LOAD_BALANCE) !!\n",
                       lv, PID, SonPID );
#           endif // LOAD_BALANCE
#           endif // DEBUG_HDF5

//          get the SonGID by "father corner = son corner -> son LB_Idx -> son GID"
//          --> didn't assume any relation between son's and father's LB_Idx
//          (although for Hilbert curve we have "SonLBIdx-SonLBIdx%8 = 8*MyLBIdx")
            SonCr    = amr->patch[0][lv][PID]->corner;
            SonLBIdx = LB_Corner2Index( SonLv, SonCr, CHECK_ON );

#           if ( defined DEBUG_HDF5  &&  LOAD_BALANCE == HILBERT )
            if ( SonLBIdx - SonLBIdx%8 != 8*amr->patch[0][lv][PID]->LB_Idx )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonCr (%d,%d,%d), incorret SonLBIdx %ld, (MyLBIdx %ld) !!\n",
                       lv, PID, SonPID, SonCr[0], SonCr[1], SonCr[2], SonLBIdx, amr->patch[0][lv][PID]->LB_Idx );
#           endif

            Mis_Matching_int( NPatchTotal[SonLv], LBIdxList_Sort[SonLv], 1, &SonLBIdx, &MatchIdx );

#           ifdef DEBUG_HDF5
            if ( MatchIdx < 0 )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d, SonLBIdx %ld, couldn't find a matching patch !!\n",
                       lv, PID, SonPID, SonLBIdx );
#           endif

            SonGID = LBIdxList_Sort_IdxTable[SonLv][MatchIdx] + GID_LvStart[SonLv];
         } // else if ( SonPID < -1 )

//       son patch is a buffer patch (SonPID >= amr->NPatchComma[SonLv][1]) --> impossible
         else // ( SonPID >= amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonPID %d is a buffer patch (NRealSonPatch %d) !!\n",
                       lv, PID, SonPID, amr->NPatchComma[SonLv][1] );

         SonList_Local[lv][PID] = SonGID;


//       4-3-5. sibling GID
         for (int s=0; s<26; s++)
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

//          no sibling (SibPID can be either -1 or SIB_OFFSET_NONPERIODIC-BoundaryDirection)
            if      ( SibPID < 0 )
               SibGID = SibPID;

//          sibling patch is a real patch
            else if ( SibPID < amr->NPatchComma[lv][1] )
               SibGID = SibPID + GID_Offset[lv];

//          sibling patch is a buffer patch (which may lie outside the simulation domain)
            else
            {
#              ifdef DEBUG_HDF5
               if ( SibPID >= amr->num[lv] )
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, SibPID %d >= total number of patches %d !!\n",
                          lv, PID, SibPID, amr->num[lv] );
#              endif

//             get the SibGID by "sibling corner -> sibling LB_Idx -> sibling GID"
               SibCr    = amr->patch[0][lv][SibPID]->corner;
               SibLBIdx = LB_Corner2Index( lv, SibCr, CHECK_OFF );   // periodicity has been assumed here

               Mis_Matching_int( NPatchTotal[lv], LBIdxList_Sort[lv], 1, &SibLBIdx, &MatchIdx );

#              ifdef DEBUG_HDF5
               if ( MatchIdx < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, SibPID %d, SibLBIdx %ld, couldn't find a matching patch !!\n",
                          lv, PID, SibPID, SibLBIdx );
#              endif

               SibGID = LBIdxList_Sort_IdxTable[lv][MatchIdx] + GID_LvStart[lv];
            } // if ( SibPID >= amr->NPatchComma[lv][1] )

            SibList_Local[lv][PID][s] = SibGID;

         } // for (int s=0; s<26; s++)


#        ifdef PARTICLE
//       4-3-6. NPar
         NParList_Local[lv][PID] = amr->patch[0][lv][PID]->NPar;
#        endif
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


// 4-4. gather data from all ranks
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_Fa  [r] = NPatchAllRank[r][lv];
         RecvCount_Son [r] = RecvCount_Fa[r];
         RecvCount_Sib [r] = RecvCount_Fa[r]*26;
         RecvCount_Cr  [r] = RecvCount_Fa[r]*3;
#        ifdef PARTICLE
         RecvCount_NPar[r] = RecvCount_Fa[r];
#        endif

         RecvDisp_Fa   [r] = ( r == 0 ) ? 0 : RecvDisp_Fa[r-1] + RecvCount_Fa[r-1];
         RecvDisp_Son  [r] = RecvDisp_Fa[r];
         RecvDisp_Sib  [r] = RecvDisp_Fa[r]*26;
         RecvDisp_Cr   [r] = RecvDisp_Fa[r]*3;
#        ifdef PARTICLE
         RecvDisp_NPar [r] = RecvDisp_Fa[r];
#        endif
      }

//    note that we collect data at one level at a time
      MPI_Gatherv( FaList_Local[lv],     amr->NPatchComma[lv][1],    MPI_INT,
                   FaList_AllLv+GID_LvStart[lv],       RecvCount_Fa,   RecvDisp_Fa,   MPI_INT, 0, MPI_COMM_WORLD );

      MPI_Gatherv( SonList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                   SonList_AllLv+GID_LvStart[lv],      RecvCount_Son,  RecvDisp_Son,  MPI_INT, 0, MPI_COMM_WORLD );

      MPI_Gatherv( SibList_Local[lv][0], amr->NPatchComma[lv][1]*26, MPI_INT,
                   (SibList_AllLv+GID_LvStart[lv])[0], RecvCount_Sib,  RecvDisp_Sib,  MPI_INT, 0, MPI_COMM_WORLD );

      MPI_Gatherv( CrList_Local[lv][0],  amr->NPatchComma[lv][1]*3,  MPI_INT,
                   (CrList_AllLv+GID_LvStart[lv])[0],  RecvCount_Cr,   RecvDisp_Cr,   MPI_INT, 0, MPI_COMM_WORLD );

#     ifdef PARTICLE
      MPI_Gatherv( NParList_Local[lv],    amr->NPatchComma[lv][1],    MPI_INT,
                   NParList_AllLv+GID_LvStart[lv],     RecvCount_NPar, RecvDisp_NPar, MPI_INT, 0, MPI_COMM_WORLD );
#     endif
   } // for (int lv=0; lv<NLEVEL; lv++)


// 4-5. dump the tree info
   if ( MPI_Rank == 0 )
   {
//    reopen file
      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

      H5_GroupID_Tree = H5Gcreate( H5_FileID, "Tree", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Tree < 0 )    Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Tree" );

//    4-5-1. LBIdx
      H5_SetDims_LBIdx = NPatchAllLv;
      H5_SpaceID_LBIdx = H5Screate_simple( 1, &H5_SetDims_LBIdx, NULL );
      H5_SetID_LBIdx   = H5Dcreate( H5_GroupID_Tree, "LBIdx", H5T_NATIVE_LONG, H5_SpaceID_LBIdx,
                                    H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_LBIdx < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "LBIdx" );

      H5_Status = H5Dwrite( H5_SetID_LBIdx, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, LBIdxList_AllLv );
      H5_Status = H5Dclose( H5_SetID_LBIdx );
      H5_Status = H5Sclose( H5_SpaceID_LBIdx );

//    4-5-2. corner
      H5_SetDims_Cr[0] = NPatchAllLv;
      H5_SetDims_Cr[1] = 3;
      H5_SpaceID_Cr    = H5Screate_simple( 2, H5_SetDims_Cr, NULL );
      H5_SetID_Cr      = H5Dcreate( H5_GroupID_Tree, "Corner", H5T_NATIVE_INT, H5_SpaceID_Cr,
                                    H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Cr < 0 )    Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Corner" );

//    attach the attribute for converting corner to physical coordinates
      H5_AttID_Cvt2Phy = H5Acreate( H5_SetID_Cr, "Cvt2Phy", H5T_NATIVE_DOUBLE, H5_SpaceID_Scalar,
                                    H5P_DEFAULT, H5P_DEFAULT );

      if ( H5_AttID_Cvt2Phy < 0 )   Aux_Error( ERROR_INFO, "failed to create the attribute \"%s\" !!\n", "Cvt2Phy" );

      H5_Status = H5Awrite( H5_AttID_Cvt2Phy, H5T_NATIVE_DOUBLE, &amr->dh[TOP_LEVEL] );
      H5_Status = H5Aclose( H5_AttID_Cvt2Phy );

      H5_Status = H5Dwrite( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Cr );
      H5_Status = H5Sclose( H5_SpaceID_Cr );

//    4-5-3. father
      H5_SetDims_Fa = NPatchAllLv;
      H5_SpaceID_Fa = H5Screate_simple( 1, &H5_SetDims_Fa, NULL );
      H5_SetID_Fa   = H5Dcreate( H5_GroupID_Tree, "Father", H5T_NATIVE_INT, H5_SpaceID_Fa,
                                 H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Fa < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Father" );

      H5_Status = H5Dwrite( H5_SetID_Fa, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, FaList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Fa );
      H5_Status = H5Sclose( H5_SpaceID_Fa );

//    4-5-4. son
      H5_SetDims_Son = NPatchAllLv;
      H5_SpaceID_Son = H5Screate_simple( 1, &H5_SetDims_Son, NULL );
      H5_SetID_Son   = H5Dcreate( H5_GroupID_Tree, "Son", H5T_NATIVE_INT, H5_SpaceID_Son,
                                  H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Son < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Son" );

      H5_Status = H5Dwrite( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Son );
      H5_Status = H5Sclose( H5_SpaceID_Son );

//    4-5-5. sibling
      H5_SetDims_Sib[0] = NPatchAllLv;
      H5_SetDims_Sib[1] = 26;
      H5_SpaceID_Sib    = H5Screate_simple( 2, H5_SetDims_Sib, NULL );
      H5_SetID_Sib      = H5Dcreate( H5_GroupID_Tree, "Sibling", H5T_NATIVE_INT, H5_SpaceID_Sib,
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_Sib < 0 )    Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Sibling" );

      H5_Status = H5Dwrite( H5_SetID_Sib, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SibList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Sib );
      H5_Status = H5Sclose( H5_SpaceID_Sib );

//    4-5-6. NPar
#     ifdef PARTICLE
      H5_SetDims_NPar = NPatchAllLv;
      H5_SpaceID_NPar = H5Screate_simple( 1, &H5_SetDims_NPar, NULL );
      H5_SetID_NPar   = H5Dcreate( H5_GroupID_Tree, "NPar", H5T_NATIVE_INT, H5_SpaceID_NPar,
                                   H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );

      if ( H5_SetID_NPar < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "NPar" );

      H5_Status = H5Dwrite( H5_SetID_NPar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, NParList_AllLv );
      H5_Status = H5Dclose( H5_SetID_NPar );
      H5_Status = H5Sclose( H5_SpaceID_NPar );
#     endif

//    close file
      H5_Status = H5Gclose( H5_GroupID_Tree );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )



// 5. output the simulation grid data (density, momentum, ... etc)
   const int FieldSizeOnePatch = sizeof(real)*CUBE(PS1);
   int  NFieldOut;
   char (*FieldName)[MAX_STRING]     = NULL;
   real (*FieldData)[PS1][PS1][PS1]  = NULL;
#  ifdef SRHD
   real (*Dens) [PS1][PS1][PS1]      = NULL;
   real (*MomX) [PS1][PS1][PS1]      = NULL;
   real (*MomY) [PS1][PS1][PS1]      = NULL;
   real (*MomZ) [PS1][PS1][PS1]      = NULL;
   real (*Engy) [PS1][PS1][PS1]      = NULL;
   real (*Temp) [PS1][PS1][PS1]      = NULL;
   real Cons[NCOMP_FLUID];
#  endif


#  ifdef MHD
   const int FCMagSizeOnePatch = sizeof(real)*PS1P1*SQR(PS1);
   char FCMagName[NCOMP_MAG][MAX_STRING];
   real (*FCMagData)[PS1P1*SQR(PS1)] = NULL;
#  endif

// 5-0. determine variable indices
   NFieldOut = NCOMP_TOTAL;

#  ifdef GRAVITY
   int PotDumpIdx = -1;
   if ( OPT__OUTPUT_POT )  PotDumpIdx = NFieldOut ++;
#  endif

#  ifdef PARTICLE
   int ParDensDumpIdx = -1;
   if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )   ParDensDumpIdx = NFieldOut ++;
#  endif

#  ifdef SRHD
   int TempDumpIdx = -1; // dump index for temperature
   TempDumpIdx = NFieldOut ++;
#  endif

#  ifdef MHD
   int CCMagDumpIdx = -1;
   if ( OPT__OUTPUT_CC_MAG )
   {
      CCMagDumpIdx = NFieldOut;
      NFieldOut   += NCOMP_MAG;
   }
#  endif


// 5-1. set the output field names
   FieldName = new char [NFieldOut][MAX_STRING];

   for (int v=0; v<NCOMP_TOTAL; v++)   sprintf( FieldName[v], FieldLabel[v] );

#  ifdef GRAVITY
   if ( OPT__OUTPUT_POT )     sprintf( FieldName[PotDumpIdx], PotLabel );
#  endif

#  ifdef PARTICLE
   if      ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_PAR_ONLY )   sprintf( FieldName[ParDensDumpIdx], "ParDens" );
   else if ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_TOTAL )      sprintf( FieldName[ParDensDumpIdx], "TotalDens" );
#  endif

#  ifdef MHD
   if ( OPT__OUTPUT_CC_MAG )
   {
      sprintf( FieldName[ CCMagDumpIdx + MAGX ], "CCMagX" );
      sprintf( FieldName[ CCMagDumpIdx + MAGY ], "CCMagY" );
      sprintf( FieldName[ CCMagDumpIdx + MAGZ ], "CCMagZ" );
   }

   for (int v=0; v<NCOMP_MAG; v++)  sprintf( FCMagName[v], MagLabel[v] );
#  endif

#  ifdef SRHD
    sprintf( FieldName[TempDumpIdx], "Temp" );
#  endif


// 5-2. initialize the "GridData" group and the datasets of all fields and magnetic field
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PS1;
   H5_SetDims_Field[2] = PS1;
   H5_SetDims_Field[3] = PS1;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
      H5_SetDims_FCMag[0] = NPatchAllLv;
      for (int t=1; t<4; t++)
      H5_SetDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_SpaceID_FCMag[v] = H5Screate_simple( 4, H5_SetDims_FCMag, NULL );
      if ( H5_SpaceID_FCMag[v] < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n",
                                                  "H5_SpaceID_FCMag", v );
   }
#  endif

   if ( MPI_Rank == 0 )
   {
//    HDF5 file must be synchronized before being written by the next rank
      SyncHDF5File( FileName );

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    create the "GridData" group
      H5_GroupID_GridData = H5Gcreate( H5_FileID, "GridData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "GridData" );

//    create the datasets of all fields
      for (int v=0; v<NFieldOut; v++)
      {
         H5_SetID_Field = H5Dcreate( H5_GroupID_GridData, FieldName[v], H5T_GAMER_REAL, H5_SpaceID_Field,
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", FieldName[v] );
         H5_Status = H5Dclose( H5_SetID_Field );
      }

//    create the datasets of all magnetic field components
#     ifdef MHD
      for (int v=0; v<NCOMP_MAG; v++)
      {
         H5_SetID_FCMag = H5Dcreate( H5_GroupID_GridData, FCMagName[v], H5T_GAMER_REAL, H5_SpaceID_FCMag[v],
                                     H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_FCMag < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", FCMagName[v] );
         H5_Status = H5Dclose( H5_SetID_FCMag );
      }
#     endif

//    close the file and group
      H5_Status = H5Gclose( H5_GroupID_GridData );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )


// 5-3. start to dump data (serial instead of parallel)
#  ifdef PARTICLE
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const bool TimingSendPar_No  = false;
   const bool PredictParPos_No  = false;   // particles synchronization is done in "Flu_CorrAfterAllSync()"
   const bool JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool SibBufPatch       = true;
   const bool FaSibBufPatch     = true;
#  else
   const bool SibBufPatch       = NULL_BOOL;
   const bool FaSibBufPatch     = NULL_BOOL;
#  endif

   int *PID0List = NULL;
#  endif // #ifdef PARTICLE

// output one level at a time so that data at the same level are consecutive on disk (even for multiple ranks)
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    5-3-0. initialize the particle density array (rho_ext) and collect particles from higher levels for outputting particle density
#     ifdef PARTICLE
      if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
      {
         Prepare_PatchData_InitParticleDensityArray( lv );

         Par_CollectParticle2OneLevel( lv, PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                       TimingSendPar_No );
      }
#     endif

      for (int TRank=0; TRank<MPI_NRank; TRank++)
      {
         if ( MPI_Rank == TRank )
         {
//          HDF5 file must be synchronized before being written by the next rank
            SyncHDF5File( FileName );

//          reopen the file and group
            H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
            if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

            H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
            if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );


//          5-3-1. dump cell-centered data
//          5-3-1-1. determine the memory space
            H5_MemDims_Field[0] = amr->NPatchComma[lv][1];
            H5_MemDims_Field[1] = PS1;
            H5_MemDims_Field[2] = PS1;
            H5_MemDims_Field[3] = PS1;

            H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
            if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_Field" );


//          5-3-1-2. determine the subset of the dataspace
            H5_Offset_Field[0] = GID_Offset[lv];
            H5_Offset_Field[1] = 0;
            H5_Offset_Field[2] = 0;
            H5_Offset_Field[3] = 0;

            H5_Count_Field [0] = amr->NPatchComma[lv][1];
            H5_Count_Field [1] = PS1;
            H5_Count_Field [2] = PS1;
            H5_Count_Field [3] = PS1;

            H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
            if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );


//          output one field at one level in one rank at a time
            FieldData = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
#           ifdef SRHD
            Dens      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
            MomX      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
            MomY      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
            MomZ      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
            Engy      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
            Temp      = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
#           endif

            for (int v=0; v<NFieldOut; v++)
            {
//             5-3-1-3. collect the target field from all patches at the current target level
//             a. gravitational potential
#              ifdef GRAVITY
               if ( v == PotDumpIdx )
               {
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     memcpy( FieldData[PID], amr->patch[ amr->PotSg[lv] ][lv][PID]->pot, FieldSizeOnePatch );
               }
               else
#              endif

//             b. particle density on grids
#              ifdef PARTICLE
               if ( v == ParDensDumpIdx )
               {
                  PID0List = new int [ amr->NPatchComma[lv][1]/8 ];

                  for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

//                we do not check minimum density here (just because it's unnecessary)
                  Prepare_PatchData( lv, Time[lv], FieldData[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List,
                                     ( OPT__OUTPUT_PAR_DENS == PAR_OUTPUT_DENS_PAR_ONLY ) ? _PAR_DENS : _TOTAL_DENS, _NONE,
                                     OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                                     MinDens_No, MinPres_No, DE_Consistency_No );

                  delete [] PID0List;
               }
               else
#              endif

//             c. cell-centered magnetic field
#              ifdef MHD
               if ( v >= CCMagDumpIdx  &&  v < CCMagDumpIdx+NCOMP_MAG )
               {
                  const int Bv = v - CCMagDumpIdx;
                  real CCMag_1Cell[NCOMP_MAG];

                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
//                   actually we only need CCMag_1Cell[Bv] here
//                   --> but the overhead of computing the other two B components is probably acceptable
                     MHD_GetCellCenteredBFieldInPatch( CCMag_1Cell, lv, PID, i, j, k, amr->MagSg[lv] );
                     FieldData[PID][k][j][i] = CCMag_1Cell[Bv];
                  }
               }
               else
#              endif

//             d. fluid variables
               {
                  for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                     memcpy( FieldData[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
#                 ifdef SRHD
                  switch (v)
                  {
                     case DENS:
                       for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                         memcpy( Dens[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
                       break;
                     case MOMX:
                       for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                         memcpy( MomX[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
                       break;
                     case MOMY:
                       for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                         memcpy( MomY[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
                       break;
                     case MOMZ:
                       for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                         memcpy( MomZ[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
                       break;
                     case ENGY:
                       for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                         memcpy( Engy[PID], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v], FieldSizeOnePatch );
                       break;
                  }
#                 endif
               }

#              ifdef SRHD
               } // for (int v=0; v<NFieldOut; v++)

//             convert conserved quantities to temperature
               for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	           for ( int i=0;i<PS1;i++  )
	           for ( int j=0;j<PS1;j++  )
	           for ( int k=0;k<PS1;k++  )
               {
	              Cons[0] = Dens[PID][i][j][k];
	              Cons[1] = MomX[PID][i][j][k];
	              Cons[2] = MomY[PID][i][j][k];
	              Cons[3] = MomZ[PID][i][j][k];
	              Cons[4] = Engy[PID][i][j][k];
#                 ifdef CHECK_FAILED_CELL_IN_FLUID
	              if(SRHD_CheckUnphysical(Cons, NULL, __FUNCTION__, __LINE__, true)) exit(EXIT_FAILURE);
#                 endif
                  Temp[PID][i][j][k] =  Hydro_Con2Temp( Cons[0], Cons[1], Cons[2], Cons[3], Cons[4], 
                                                        NULL, NULL_BOOL, NULL_REAL, NULL_REAL, NULL,
                                                        EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
               }

//             copy conserved data and temperature into FieldData
               for (int v=0; v<NFieldOut; v++)
	           {
	              switch (v)
	               {
	                  case DENS:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], Dens[PID], FieldSizeOnePatch );
	                    break;
	                  case MOMX:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], MomX[PID], FieldSizeOnePatch );
	                    break;
	                  case MOMY:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], MomY[PID], FieldSizeOnePatch );
	                    break;
	                  case MOMZ:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], MomZ[PID], FieldSizeOnePatch );
	                    break;
	                  case ENGY:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], Engy[PID], FieldSizeOnePatch );
	                    break;
	                  case 5:
	                     for ( int PID=0;PID < amr->NPatchComma[lv][1];PID++)
	                           memcpy( FieldData[PID], Temp[PID], FieldSizeOnePatch );
	                    break;
#                     ifdef GRAVITY
                      case 6:
                         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                            memcpy( FieldData[PID], amr->patch[ amr->PotSg[lv] ][lv][PID]->pot, FieldSizeOnePatch );
                        break;
#                     endif
	                  default:
                      Aux_Error(ERROR_INFO, "Out of the range! NFieldOut=%d, v=%d\n", NFieldOut, v);
	                    break;
                   }
#              endif

//             5-3-1-4. write data to disk
               H5_SetID_Field = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );

               H5_Status = H5Dwrite( H5_SetID_Field, H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT, FieldData );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to write a field (lv %d, v %d) !!\n", lv, v );

               H5_Status = H5Dclose( H5_SetID_Field );
            } // for (int v=0; v<NFieldOut; v++)


//          5-3-1-5.free resource before dumping magnetic field to save memory
            delete [] FieldData;
#           ifdef SRHD
            delete [] Dens;
            delete [] MomX;
            delete [] MomY;
            delete [] MomZ;
            delete [] Engy;
            delete [] Temp;
#           endif

            H5_Status = H5Sclose( H5_MemID_Field );

//          free memory used for outputting particle density
#           ifdef PARTICLE
            if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE )
            {
               Prepare_PatchData_FreeParticleDensityArray( lv );

               Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );
            }
#           endif


//          5-3-2. dump magnetic field
#           ifdef MHD
//          5-3-2-0. allocate memory
//                   --> output one B component at one level in one rank at a time
            FCMagData = new real [ amr->NPatchComma[lv][1] ][ PS1P1*SQR(PS1) ];

            for (int v=0; v<NCOMP_MAG; v++)
            {
//             5-3-2-1. determine the memory space
               H5_MemDims_FCMag[0] = amr->NPatchComma[lv][1];
               for (int t=1; t<4; t++)
               H5_MemDims_FCMag[t] = ( 3-t == v ) ? PS1P1 : PS1;

               H5_MemID_FCMag = H5Screate_simple( 4, H5_MemDims_FCMag, NULL );
               if ( H5_MemID_FCMag < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_FCMag" );


//             5-3-2-2. determine the subset of the dataspace
               H5_Offset_FCMag[0] = GID_Offset[lv];
               H5_Offset_FCMag[1] = 0;
               H5_Offset_FCMag[2] = 0;
               H5_Offset_FCMag[3] = 0;

               H5_Count_FCMag [0] = amr->NPatchComma[lv][1];
               for (int t=1; t<4; t++)
               H5_Count_FCMag [t] = ( 3-t == v ) ? PS1P1 : PS1;

               H5_Status = H5Sselect_hyperslab( H5_SpaceID_FCMag[v], H5S_SELECT_SET, H5_Offset_FCMag, NULL, H5_Count_FCMag, NULL );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the magnetic field !!\n" );


//             5-3-2-3. collect the target B component from all patches at the current target level
               for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
                  memcpy( FCMagData[PID], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v], FCMagSizeOnePatch );


//             5-3-2-4. write data to disk
               H5_SetID_FCMag = H5Dopen( H5_GroupID_GridData, FCMagName[v], H5P_DEFAULT );

               H5_Status = H5Dwrite( H5_SetID_FCMag, H5T_GAMER_REAL, H5_MemID_FCMag, H5_SpaceID_FCMag[v], H5P_DEFAULT, FCMagData );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to write magnetic field (lv %d, v %d) !!\n", lv, v );

               H5_Status = H5Dclose( H5_SetID_FCMag );
               H5_Status = H5Sclose( H5_MemID_FCMag );
            } // for (int v=0; v<NCOMP_MAG; v++)


//          5-3-2-5.free resource
            delete [] FCMagData;
#           endif // #ifdef MHD

            H5_Status = H5Gclose( H5_GroupID_GridData );
            H5_Status = H5Fclose( H5_FileID );
         } // if ( MPI_Rank == TRank )

         MPI_Barrier( MPI_COMM_WORLD );

      } // for (int TRank=0; TRank<MPI_NRank; TRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)

   H5_Status = H5Sclose( H5_SpaceID_Field );
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   H5_Status = H5Sclose( H5_SpaceID_FCMag[v] );
#  endif



// 6. output particles
#  ifdef PARTICLE
//###ISSUE: currently we output all particles at the same level at once (although one attribute at a time),
//          which may introduce a large memory overhead
//          --> solution: we can output a fixed number of particles at a time (see Output_DumpData_Total.cpp)
   long (*NParLv_EachRank)[NLEVEL] = new long [MPI_NRank][NLEVEL];   // number of particles at each level in each rank
   real (*ParBuf1v1Lv)             = NULL;   // buffer storing the data of one particle attribute at one level

   long  GParID_Offset[NLEVEL];  // GParID = global particle index (==> unique for each particle)
   long  NParLv_AllRank[NLEVEL];
   long  MaxNPar1Lv, NParInBuf, ParID;


// 6-1. initialize variables
// 6-1-1. allocate I/O buffer for storing particle data
   MaxNPar1Lv = 0;
   for (int lv=0; lv<NLEVEL; lv++)  MaxNPar1Lv = MAX( MaxNPar1Lv, amr->Par->NPar_Lv[lv] );

   ParBuf1v1Lv = new real [MaxNPar1Lv];

// 6-1-2. get the starting global particle index (i.e., GParID_Offset[NLEVEL]) for particles at each level in this rank
   MPI_Allgather( amr->Par->NPar_Lv, NLEVEL, MPI_LONG, NParLv_EachRank[0], NLEVEL, MPI_LONG, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      NParLv_AllRank[lv] = 0;
      for (int r=0; r<MPI_NRank; r++)     NParLv_AllRank[lv] += NParLv_EachRank[r][lv];

      GParID_Offset[lv] = 0;
      for (int FaLv=0; FaLv<lv; FaLv++)   GParID_Offset[lv] += NParLv_AllRank[FaLv];

      for (int r=0; r<MPI_Rank; r++)      GParID_Offset[lv] += NParLv_EachRank[r][lv];
   }


// 6-2. initialize the "Particle" group and the datasets of all particle attributes
   H5_SetDims_ParData[0] = amr->Par->NPar_Active_AllRank;
   H5_SpaceID_ParData    = H5Screate_simple( 1, H5_SetDims_ParData, NULL );
   if ( H5_SpaceID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_ParData" );

   if ( MPI_Rank == 0 )
   {
//    HDF5 file must be synchronized before being written by the next rank
      SyncHDF5File( FileName );

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    create the "Particle" group
      H5_GroupID_Particle = H5Gcreate( H5_FileID, "Particle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Particle" );

//    create the datasets of all particle attributes
      for (int v=0; v<PAR_NATT_STORED; v++)
      {
         H5_SetID_ParData = H5Dcreate( H5_GroupID_Particle, ParAttLabel[v], H5T_GAMER_REAL, H5_SpaceID_ParData,
                                       H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
         if ( H5_SetID_ParData < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", ParAttLabel[v] );
         H5_Status = H5Dclose( H5_SetID_ParData );
      }

//    close the file and group
      H5_Status = H5Gclose( H5_GroupID_Particle );
      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )


// 6-3. start to dump particle data (one level, one rank, and one attribute at a time)
//      --> note that particles must be outputted in the same order as their associated patches
   for (int lv=0; lv<NLEVEL; lv++)
   for (int TRank=0; TRank<MPI_NRank; TRank++)
   {
      if ( MPI_Rank == TRank )
      {
//       HDF5 file must be synchronized before being written by the next rank
         SyncHDF5File( FileName );

//       reopen the file and group
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDWR, H5P_DEFAULT );
         if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_Particle = H5Gopen( H5_FileID, "Particle", H5P_DEFAULT );
         if ( H5_GroupID_Particle < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Particle" );


//       6-3-1. determine the memory space
         H5_MemDims_ParData[0] = amr->Par->NPar_Lv[lv];
         H5_MemID_ParData      = H5Screate_simple( 1, H5_MemDims_ParData, NULL );
         if ( H5_MemID_ParData < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_ParData" );


//       6-3-2. determine the subset of the dataspace
         H5_Offset_ParData[0] = GParID_Offset[lv];
         H5_Count_ParData [0] = amr->Par->NPar_Lv[lv];

         H5_Status = H5Sselect_hyperslab( H5_SpaceID_ParData, H5S_SELECT_SET, H5_Offset_ParData, NULL, H5_Count_ParData, NULL );
         if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the particle data !!\n" );


//       output one particle attribute at one level in one rank at a time
//       --> skip the last PAR_NATT_UNSTORED attributes since we do not want to store them on disk
         for (int v=0; v<PAR_NATT_STORED; v++)
         {
//          6-3-3. collect particle data from all patches at the current target level
            NParInBuf = 0;

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
            {
               ParID = amr->patch[0][lv][PID]->ParList[p];

#              ifdef DEBUG_PARTICLE
               if ( NParInBuf >= amr->Par->NPar_Lv[lv] )
                  Aux_Error( ERROR_INFO, "lv %d, NParInBuf (%ld) >= NPar_Lv (%ld) !!\n", lv, NParInBuf, amr->Par->NPar_Lv[lv] );
#              endif

               ParBuf1v1Lv[ NParInBuf ++ ] = amr->Par->Attribute[v][ParID];
            }


//          6-3-4. write data to disk
            H5_SetID_ParData = H5Dopen( H5_GroupID_Particle, ParAttLabel[v], H5P_DEFAULT );

            H5_Status = H5Dwrite( H5_SetID_ParData, H5T_GAMER_REAL, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT, ParBuf1v1Lv );
            if ( H5_Status < 0 )
               Aux_Error( ERROR_INFO, "failed to write a particle attribute (lv %d, v %d) !!\n", lv, v );

            H5_Status = H5Dclose( H5_SetID_ParData );
         } // for (int v=0; v<PAR_NATT_STORED; v++)

//       free resource
         H5_Status = H5Sclose( H5_MemID_ParData );
         H5_Status = H5Gclose( H5_GroupID_Particle );
         H5_Status = H5Fclose( H5_FileID );
      } // if ( MPI_Rank == TRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TRank=0; TRank<MPI_NRank; TRank++) ... for (int lv=0; lv<NLEVEL; lv++)

   H5_Status = H5Sclose( H5_SpaceID_ParData );

   delete [] ParBuf1v1Lv;
   delete [] NParLv_EachRank;
#  endif // #ifdef PARTICLE



// 7. check
#  ifdef DEBUG_HDF5
   if ( MPI_Rank == 0 )
   {
      const int MirrorSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };

      H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
      if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName );

//    7-1. validate the father-son relation
//    7-1-1. load data
      char SetName[MAX_STRING];
      sprintf( SetName, "Tree/Father" );
      H5_SetID_Fa = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Fa < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Fa, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, FaList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Fa );

      sprintf( SetName, "Tree/Son" );
      H5_SetID_Son = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Son < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Son );

      sprintf( SetName, "Tree/Sibling" );
      H5_SetID_Sib = H5Dopen( H5_FileID, SetName, H5P_DEFAULT );

      if ( H5_SetID_Sib < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", SetName );

      H5_Status = H5Dread( H5_SetID_Sib, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SibList_AllLv );
      H5_Status = H5Dclose( H5_SetID_Sib );


      for (int lv=0; lv<NLEVEL; lv++)
      for (int GID=GID_LvStart[lv]; GID<GID_LvStart[lv]+NPatchTotal[lv]; GID++)
      {
//       7-1-2. root patches have no father
         if ( lv == 0 )
         if ( FaList_AllLv[GID] != -1 )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d != -1 !!\n", lv, GID, FaList_AllLv[GID] );

//       7-1-3. all patches at refinement levels have fathers
         if ( lv > 0 )
         if ( FaList_AllLv[GID] < 0  ||  FaList_AllLv[GID] >= GID_LvStart[lv] )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d < 0 (or > max = %d) !!\n",
                       lv, GID, FaList_AllLv[GID], GID_LvStart[lv]-1 );

//       7-1-4. father->son == itself
         if ( lv > 0 )
         if ( SonList_AllLv[ FaList_AllLv[GID] ] + GID%8 != GID )
            Aux_Error( ERROR_INFO, "Lv %d, GID %d, FaGID %d, FaGID->Son %d ==> inconsistent !!\n",
                       lv, GID, FaList_AllLv[GID], SonList_AllLv[ FaList_AllLv[GID] ] );

//       7-1-5. son->father == itself
         SonGID = SonList_AllLv[GID];
         if ( SonGID != -1 )
         {
            if ( lv >= MAX_LEVEL )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d != -1 !!\n", lv, GID, SonGID );

            if ( SonGID < -1 )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d < -1 !!\n", lv, GID, SonGID );

            if ( lv < NLEVEL-1  &&  SonGID >= GID_LvStart[lv+1]+NPatchTotal[lv+1] )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d > max %d !!\n", lv, GID, SonGID,
                          GID_LvStart[lv+1]+NPatchTotal[lv+1]-1 );

            for (int LocalID=0; LocalID<8; LocalID++)
            if ( FaList_AllLv[SonGID+LocalID] != GID )
               Aux_Error( ERROR_INFO, "Lv %d, GID %d, SonGID %d, SonGID->Father %d ==> inconsistent !!\n",
                          lv, GID, SonGID+LocalID, FaList_AllLv[SonGID+LocalID] );
         }

//       7-1-6. sibling->sibling_mirror = itself
         for (int s=0; s<26; s++)
         {
            SibGID = SibList_AllLv[GID][s];

            if ( SibGID >= 0 )
            {
               if ( SibGID < GID_LvStart[lv]  ||  SibGID >= GID_LvStart[lv]+NPatchTotal[lv] )
                  Aux_Error( ERROR_INFO, "Lv %d, GID %d, sib %d, SibGID %d lies outside the correct range (%d <= SibGID < %d) !!\n",
                             lv, GID, s, SibGID, GID_LvStart[lv], GID_LvStart[lv]+NPatchTotal[lv] );

               if ( SibList_AllLv[SibGID][ MirrorSib[s] ] != GID )
                  Aux_Error( ERROR_INFO, "Lv %d, GID %d, sib %d, SibGID %d != SibGID->sibling %d !!\n",
                             lv, GID, s, SibGID, SibList_AllLv[SibGID][ MirrorSib[s] ] );
            }
         }
      } // for (int lv=0; lv<NLEVEL; lv++)

      H5_Status = H5Fclose( H5_FileID );
   } // if ( MPI_Rank == 0 )
#  endif // #ifdef DEBUG_HDF5



// 8. close all HDF5 objects and free memory
   H5_Status = H5Tclose( H5_TypeID_Com_KeyInfo );
   H5_Status = H5Tclose( H5_TypeID_Com_Makefile );
   H5_Status = H5Tclose( H5_TypeID_Com_SymConst );
   H5_Status = H5Tclose( H5_TypeID_Com_InputPara );
   H5_Status = H5Sclose( H5_SpaceID_Scalar );
   H5_Status = H5Pclose( H5_DataCreatePropList );

   delete [] NPatchAllRank;
   delete [] FieldName;

   if ( MPI_Rank == 0 )
   {
      delete [] LBIdxList_AllLv;
      delete []    CrList_AllLv;
      delete []    FaList_AllLv;
      delete []   SonList_AllLv;
      delete []   SibList_AllLv;
#     ifdef PARTICLE
      delete []  NParList_AllLv;
#     endif
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete [] LBIdxList_Local[lv];
      delete []    CrList_Local[lv];
      delete []    FaList_Local[lv];
      delete []   SonList_Local[lv];
      delete []   SibList_Local[lv];
#     ifdef PARTICLE
      delete []  NParList_Local[lv];
#     endif

      delete [] LBIdxList_Sort[lv];
      delete [] LBIdxList_Sort_IdxTable[lv];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Total_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_KeyInfo
// Description :  Fill in the KeyInfo_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  KeyInfo : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_KeyInfo( KeyInfo_t &KeyInfo )
{

   const time_t CalTime = time( NULL );   // calendar time

   KeyInfo.FormatVersion        = 2429;
   KeyInfo.Model                = MODEL;
   KeyInfo.NLevel               = NLEVEL;
   KeyInfo.NCompFluid           = NCOMP_FLUID;
   KeyInfo.NCompPassive         = NCOMP_PASSIVE;
   KeyInfo.PatchSize            = PS1;
   KeyInfo.DumpID               = DumpID;
   KeyInfo.Step                 = Step;
#  ifdef GRAVITY
   KeyInfo.AveDens_Init         = AveDensity_Init;
   KeyInfo.Gravity              = 1;
#  else
   KeyInfo.Gravity              = 0;
#  endif
#  ifdef PARTICLE
   KeyInfo.Particle             = 1;
#  else
   KeyInfo.Particle             = 0;
#  endif
#  ifdef FLOAT8
   KeyInfo.Float8               = 1;
#  else
   KeyInfo.Float8               = 0;
#  endif
#  ifdef PARTICLE
   KeyInfo.Par_NPar             = amr->Par->NPar_Active_AllRank;
   KeyInfo.Par_NAttStored       = PAR_NATT_STORED;
#  endif
#  if ( MODEL == HYDRO )
#  ifdef MHD
   KeyInfo.Magnetohydrodynamics = 1;
#  else
   KeyInfo.Magnetohydrodynamics = 0;
#  endif
#  ifdef COSMIC_RAY
   KeyInfo.CosmicRay            = 1;
#  else
   KeyInfo.CosmicRay            = 0;
#  endif
#  endif // #if ( MODEL == HYDRO )

   for (int d=0; d<3; d++)
   {
      KeyInfo.NX0     [d] = NX0_TOT      [d];
      KeyInfo.BoxScale[d] = amr->BoxScale[d];
      KeyInfo.BoxSize [d] = amr->BoxSize [d];
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      KeyInfo.Time          [lv] = Time          [lv];
      KeyInfo.CellSize      [lv] = amr->dh       [lv];
      KeyInfo.CellScale     [lv] = amr->scale    [lv];
      KeyInfo.NPatch        [lv] = NPatchTotal   [lv];
      KeyInfo.AdvanceCounter[lv] = AdvanceCounter[lv];
      KeyInfo.dTime_AllLv   [lv] = dTime_AllLv   [lv];
   }

   KeyInfo.CodeVersion  = (char*)VERSION;
   KeyInfo.DumpWallTime = ctime( &CalTime );
   KeyInfo.DumpWallTime[ strlen(KeyInfo.DumpWallTime)-1 ] = '\0';  // remove the last character '\n'

} // FUNCTION : FillIn_KeyInfo



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_Makefile
// Description :  Fill in the Makefile_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  Makefile : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_Makefile( Makefile_t &Makefile )
{

// model-independent options
   Makefile.Model                  = MODEL;

#  ifdef GRAVITY
   Makefile.Gravity                = 1;
#  else
   Makefile.Gravity                = 0;
#  endif

#  ifdef COMOVING
   Makefile.Comoving               = 1;
#  else
   Makefile.Comoving               = 0;
#  endif

#  ifdef PARTICLE
   Makefile.Particle               = 1;
#  else
   Makefile.Particle               = 0;
#  endif


#  ifdef GPU
   Makefile.UseGPU                 = 1;
#  else
   Makefile.UseGPU                 = 0;
#  endif

#  ifdef GAMER_DEBUG
   Makefile.GAMER_Debug            = 1;
#  else
   Makefile.GAMER_Debug            = 0;
#  endif

#  ifdef BITWISE_REPRODUCIBILITY
   Makefile.BitwiseReproducibility = 1;
#  else
   Makefile.BitwiseReproducibility = 0;
#  endif

#  ifdef TIMING
   Makefile.Timing                 = 1;
#  else
   Makefile.Timing                 = 0;
#  endif

#  ifdef TIMING_SOLVER
   Makefile.TimingSolver           = 1;
#  else
   Makefile.TimingSolver           = 0;
#  endif

#  ifdef FLOAT8
   Makefile.Float8                 = 1;
#  else
   Makefile.Float8                 = 0;
#  endif

#  ifdef SERIAL
   Makefile.Serial                 = 1;
#  else
   Makefile.Serial                 = 0;
#  endif

#  ifdef LOAD_BALANCE
   Makefile.LoadBalance            = LOAD_BALANCE;
#  else
   Makefile.LoadBalance            = 0;
#  endif

#  ifdef OVERLAP_MPI
   Makefile.OverlapMPI             = 1;
#  else
   Makefile.OverlapMPI             = 0;
#  endif

#  ifdef OPENMP
   Makefile.OpenMP                 = 1;
#  else
   Makefile.OpenMP                 = 0;
#  endif

#  ifdef GPU
   Makefile.GPU_Arch               = GPU_ARCH;
#  else
   Makefile.GPU_Arch               = NULL_INT;
#  endif

#  ifdef LAOHU
   Makefile.Laohu                  = 1;
#  else
   Makefile.Laohu                  = 0;
#  endif

#  ifdef SUPPORT_HDF5
   Makefile.SupportHDF5            = 1;
#  else
   Makefile.SupportHDF5            = 0;
#  endif

#  ifdef SUPPORT_GSL
   Makefile.SupportGSL             = 1;
#  else
   Makefile.SupportGSL             = 0;
#  endif

#  ifdef SUPPORT_GRACKLE
   Makefile.SupportGrackle         = 1;
#  else
   Makefile.SupportGrackle         = 0;
#  endif

   Makefile.RandomNumber           = RANDOM_NUMBER;

   Makefile.NLevel                 = NLEVEL;
   Makefile.MaxPatch               = MAX_PATCH;


// model-dependent options
#  ifdef GRAVITY
   Makefile.PotScheme              = POT_SCHEME;

#  ifdef STORE_POT_GHOST
   Makefile.StorePotGhost          = 1;
#  else
   Makefile.StorePotGhost          = 0;
#  endif

#  ifdef UNSPLIT_GRAVITY
   Makefile.UnsplitGravity         = 1;
#  else
   Makefile.UnsplitGravity         = 0;
#  endif

#  endif // #ifdef GRAVITY

#  if   ( MODEL == HYDRO )
   Makefile.FluScheme              = FLU_SCHEME;

#  ifdef LR_SCHEME
   Makefile.LRScheme               = LR_SCHEME;
#  endif

#  ifdef RSOLVER
   Makefile.RSolver                = RSOLVER;
#  endif

#  ifdef DUAL_ENERGY
   Makefile.DualEnergy             = DUAL_ENERGY;
#  else
   Makefile.DualEnergy             = 0;
#  endif

#  ifdef MHD
   Makefile.Magnetohydrodynamics   = 1;
#  else
   Makefile.Magnetohydrodynamics   = 0;
#  endif

#  ifdef COSMIC_RAY
   Makefile.CosmicRay              = 1;
#  else
   Makefile.CosmicRay              = 0;
#  endif

   Makefile.EoS                    = EOS;

#  ifdef SRHD
   Makefile.SRHydrodynamics        = 1;
   Makefile.SRHydroReducedEnergy   = REDUCED_ENERGY;
   Makefile.SRHydroFourVelocity    = FOUR_VELOCITY;
#  else
   Makefile.SRHydrodynamics        = 0;
#  endif

#  ifdef BAROTROPIC_EOS
   Makefile.BarotropicEoS          = 1;
#  else
   Makefile.BarotropicEoS          = 0;
#  endif


#  elif ( MODEL == ELBDM )
#  ifdef CONSERVE_MASS
   Makefile.ConserveMass           = 1;
#  else
   Makefile.ConserveMass           = 0;
#  endif

#  ifdef LAPLACIAN_4TH
   Makefile.Laplacian4th           = 1;
#  else
   Makefile.Laplacian4th           = 0;
#  endif

#  ifdef QUARTIC_SELF_INTERACTION
   Makefile.SelfInteraction4       = 1;
#  else
   Makefile.SelfInteraction4       = 0;
#  endif


#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
#  ifdef STORE_PAR_ACC
   Makefile.StoreParAcc            = 1;
#  else
   Makefile.StoreParAcc            = 0;
#  endif

#  ifdef STAR_FORMATION
   Makefile.StarFormation          = 1;
#  else
   Makefile.StarFormation          = 0;
#  endif

   Makefile.Par_NAttUser           = PAR_NATT_USER;
#  endif

} // FUNCTION : FillIn_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_SymConst
// Description :  Fill in the SymConst_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  SymConst : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_SymConst( SymConst_t &SymConst )
{

// model-independent variables
   SymConst.NCompFluid           = NCOMP_FLUID;
   SymConst.NCompPassive         = NCOMP_PASSIVE;
   SymConst.PatchSize            = PS1;
   SymConst.Flu_NIn              = FLU_NIN;
   SymConst.Flu_NOut             = FLU_NOUT;
   SymConst.Flu_NIn_T            = FLU_NIN_T;
   SymConst.Flu_NIn_S            = FLU_NIN_S;
   SymConst.Flu_NOut_S           = FLU_NOUT_S;
   SymConst.NFluxFluid           = NFLUX_FLUID;
   SymConst.NFluxPassive         = NFLUX_PASSIVE;
   SymConst.Flu_GhostSize        = FLU_GHOST_SIZE;
   SymConst.Flu_Nxt              = FLU_NXT;
#  ifdef DEBUG_HDF5
   SymConst.Debug_HDF5           = 1;
#  else
   SymConst.Debug_HDF5           = 0;
#  endif
   SymConst.SibOffsetNonperiodic = SIB_OFFSET_NONPERIODIC;
#  ifdef LOAD_BALANCE
   SymConst.SonOffsetLB          = SON_OFFSET_LB;
#  endif

   SymConst.TinyNumber           = TINY_NUMBER;
   SymConst.HugeNumber           = HUGE_NUMBER;
   SymConst.MaxError             = MAX_ERROR;


// model-dependent variables
#  ifdef GRAVITY
   SymConst.Gra_NIn              = GRA_NIN;
   SymConst.Pot_GhostSize        = POT_GHOST_SIZE;
   SymConst.Gra_GhostSize        = GRA_GHOST_SIZE;
   SymConst.Rho_GhostSize        = RHO_GHOST_SIZE;
   SymConst.Pot_Nxt              = POT_NXT;
   SymConst.Gra_Nxt              = GRA_NXT;
   SymConst.Rho_Nxt              = RHO_NXT;

#  ifdef UNSPLIT_GRAVITY
   SymConst.USG_GhostSizeF       = USG_GHOST_SIZE_F;
   SymConst.USG_GhostSizeG       = USG_GHOST_SIZE_G;
   SymConst.USG_NxtF             = USG_NXT_F;
   SymConst.USG_NxtG             = USG_NXT_G;
#  endif

   SymConst.ExtPot_BlockSize     = EXTPOT_BLOCK_SIZE;
   SymConst.Gra_BlockSize        = GRA_BLOCK_SIZE;
   SymConst.ExtPotNAuxMax        = EXT_POT_NAUX_MAX;
   SymConst.ExtAccNAuxMax        = EXT_ACC_NAUX_MAX;

#  if   ( POT_SCHEME == SOR )
   SymConst.Pot_BlockSize_z      = POT_BLOCK_SIZE_Z;
#  ifdef USE_PSOLVER_10TO14
   SymConst.UsePSolver_10to14    = 1;
#  else
   SymConst.UsePSolver_10to14    = 0;
#  endif
#  ifdef SOR_RHO_SHARED
   SymConst.SOR_RhoShared        = 1;
#  else
   SymConst.SOR_RhoShared        = 0;
#  endif
#  ifdef SOR_CPOT_SHARED
   SymConst.SOR_CPotShared       = 1;
#  else
   SymConst.SOR_CPotShared       = 0;
#  endif
#  ifdef SOR_USE_SHUFFLE
   SymConst.SOR_UseShuffle       = 1;
#  else
   SymConst.SOR_UseShuffle       = 0;
#  endif
#  ifdef SOR_USE_PADDING
   SymConst.SOR_UsePadding       = 1;
#  else
   SymConst.SOR_UsePadding       = 0;
#  endif
   SymConst.SOR_ModReduction     = SOR_MOD_REDUCTION;

#  elif ( POT_SCHEME == MG  )
   SymConst.Pot_BlockSize_x      = POT_BLOCK_SIZE_X;
#  endif // POT_SCHEME
#  endif // #ifdef GRAVITY


#  ifdef PARTICLE
   SymConst.Par_NAttStored       = PAR_NATT_STORED;
   SymConst.RhoExt_GhostSize     = RHOEXT_GHOST_SIZE;

#  ifdef DEBUG_PARTICLE
   SymConst.Debug_Particle       = 1;
#  else
   SymConst.Debug_Particle       = 0;
#  endif

   SymConst.ParList_GrowthFactor = PARLIST_GROWTH_FACTOR;
   SymConst.ParList_ReduceFactor = PARLIST_REDUCE_FACTOR;
#  endif


#  ifdef BIT_REP_FLUX
   SymConst.BitRep_Flux          = 1;
#  else
   SymConst.BitRep_Flux          = 0;
#  endif

#  ifdef MHD
#  ifdef BIT_REP_ELECTRIC
   SymConst.BitRep_Electric      = 1;
#  else
   SymConst.BitRep_Electric      = 0;
#  endif
#  endif // #ifdef MHD


#  if   ( MODEL == HYDRO )
   SymConst.Flu_BlockSize_x      = FLU_BLOCK_SIZE_X;
   SymConst.Flu_BlockSize_y      = FLU_BLOCK_SIZE_Y;
#  ifdef CHECK_NEGATIVE_IN_FLUID
   SymConst.CheckNegativeInFluid = 1;
#  else
   SymConst.CheckNegativeInFluid = 0;
#  endif
#  ifdef CHAR_RECONSTRUCTION
   SymConst.CharReconstruction   = 1;
#  else
   SymConst.CharReconstruction   = 0;
#  endif
#  ifdef LR_EINT
   SymConst.LR_Eint              = 1;
#  else
   SymConst.LR_Eint              = 0;
#  endif
#  ifdef CHECK_INTERMEDIATE
   SymConst.CheckIntermediate    = CHECK_INTERMEDIATE;
#  else
   SymConst.CheckIntermediate    = 0;
#  endif
#  ifdef HLL_NO_REF_STATE
   SymConst.HLL_NoRefState       = 1;
#  else
   SymConst.HLL_NoRefState       = 0;
#  endif
#  ifdef HLL_INCLUDE_ALL_WAVES
   SymConst.HLL_IncludeAllWaves  = 1;
#  else
   SymConst.HLL_IncludeAllWaves  = 0;
#  endif
   SymConst.HLLC_WaveSpeed       = HLLC_WAVESPEED;
   SymConst.HLLE_WaveSpeed       = HLLE_WAVESPEED;
#  ifdef MHD
   SymConst.HLLD_WaveSpeed       = HLLD_WAVESPEED;
#  endif
#  ifdef N_FC_VAR
   SymConst.N_FC_Var             = N_FC_VAR;
#  endif
#  ifdef N_SLOPE_PPM
   SymConst.N_Slope_PPM          = N_SLOPE_PPM;
#  endif
#  ifdef MHD
#  ifdef EULERY
   SymConst.EulerY               = 1;
#  else
   SymConst.EulerY               = 0;
#  endif
#  endif // MHD
   SymConst.EoSNAuxMax           = EOS_NAUX_MAX;
   SymConst.EoSNTableMax         = EOS_NTABLE_MAX;

#  elif  ( MODEL == ELBDM )
   SymConst.Flu_BlockSize_x      = FLU_BLOCK_SIZE_X;
   SymConst.Flu_BlockSize_y      = FLU_BLOCK_SIZE_Y;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   SymConst.dt_Flu_BlockSize     = DT_FLU_BLOCK_SIZE;
#  ifdef DT_FLU_USE_SHUFFLE
   SymConst.dt_Flu_UseShuffle    = 1;
#  else
   SymConst.dt_Flu_UseShuffle    = 0;
#  endif
#  ifdef GRAVITY
   SymConst.dt_Gra_BlockSize     = DT_GRA_BLOCK_SIZE;
#  ifdef DT_GRA_USE_SHUFFLE
   SymConst.dt_Gra_UseShuffle    = 1;
#  else
   SymConst.dt_Gra_UseShuffle    = 0;
#  endif
#  endif

   SymConst.Src_BlockSize        = SRC_BLOCK_SIZE;
   SymConst.Src_GhostSize        = SRC_GHOST_SIZE;
   SymConst.Src_NAuxDlep         = SRC_NAUX_DLEP;
   SymConst.Src_DlepProfNVar     = SRC_DLEP_PROF_NVAR;
   SymConst.Src_DlepProfNBinMax  = SRC_DLEP_PROF_NBINMAX;
   SymConst.Src_NAuxUser         = SRC_NAUX_USER;

} // FUNCTION : FillIn_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  FillIn_InputPara
// Description :  Fill in the InputPara_t structure
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. Call-by-reference
//
// Parameter   :  InputPara : Pointer to be filled in
//-------------------------------------------------------------------------------------------------------
void FillIn_InputPara( InputPara_t &InputPara )
{

// simulation scale
   InputPara.BoxSize                 = BOX_SIZE;
   for (int d=0; d<3; d++)
   InputPara.NX0_Tot[d]              = NX0_TOT[d];
   InputPara.MPI_NRank               = MPI_NRank;
   for (int d=0; d<3; d++)
   InputPara.MPI_NRank_X[d]          = MPI_NRank_X[d];
   InputPara.OMP_NThread             = OMP_NTHREAD;
   InputPara.EndT                    = END_T;
   InputPara.EndStep                 = END_STEP;

// test problems
   InputPara.TestProb_ID             = TESTPROB_ID;

// code units
   InputPara.Opt__Unit               = OPT__UNIT;
   InputPara.Unit_L                  = UNIT_L;
   InputPara.Unit_M                  = UNIT_M;
   InputPara.Unit_T                  = UNIT_T;
   InputPara.Unit_V                  = UNIT_V;
   InputPara.Unit_D                  = UNIT_D;
   InputPara.Unit_E                  = UNIT_E;
   InputPara.Unit_P                  = UNIT_P;
#  ifdef MHD
   InputPara.Unit_B                  = UNIT_B;
#  endif

// boundary condition
   for (int t=0; t<6; t++)
   InputPara.Opt__BC_Flu[t]          = OPT__BC_FLU[t];
#  ifdef GRAVITY
   InputPara.Opt__BC_Pot             = OPT__BC_POT;
   InputPara.GFunc_Coeff0            = GFUNC_COEFF0;
#  endif

// particle
#  ifdef PARTICLE
   InputPara.Par_Init                = amr->Par->Init;
   InputPara.Par_ICFormat            = amr->Par->ParICFormat;
   InputPara.Par_ICMass              = amr->Par->ParICMass;
   InputPara.Par_Interp              = amr->Par->Interp;
   InputPara.Par_Integ               = amr->Par->Integ;
   InputPara.Par_ImproveAcc          = amr->Par->ImproveAcc;
   InputPara.Par_PredictPos          = amr->Par->PredictPos;
   InputPara.Par_RemoveCell          = amr->Par->RemoveCell;
   InputPara.Par_GhostSize           = amr->Par->GhostSize;
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   InputPara.ParAttLabel[v]          = ParAttLabel[v];
#  endif

// cosmology
#  ifdef COMOVING
   InputPara.A_Init                  = A_INIT;
   InputPara.OmegaM0                 = OMEGA_M0;
   InputPara.Hubble0                 = HUBBLE0;
#  endif

// time-step determination
   InputPara.Dt__Max                 = DT__MAX;
   InputPara.Dt__Fluid               = DT__FLUID;
   InputPara.Dt__FluidInit           = DT__FLUID_INIT;
#  ifdef GRAVITY
   InputPara.Dt__Gravity             = DT__GRAVITY;
#  endif
#  if ( MODEL == ELBDM )
   InputPara.Dt__Phase               = DT__PHASE;
#  endif
#  ifdef PARTICLE
   InputPara.Dt__ParVel              = DT__PARVEL;
   InputPara.Dt__ParVelMax           = DT__PARVEL_MAX;
   InputPara.Dt__ParAcc              = DT__PARACC;
#  endif
#  ifdef COMOVING
   InputPara.Dt__MaxDeltaA           = DT__MAX_DELTA_A;
#  endif
   InputPara.Dt__SyncParentLv        = DT__SYNC_PARENT_LV;
   InputPara.Dt__SyncChildrenLv      = DT__SYNC_CHILDREN_LV;
   InputPara.Opt__DtUser             = OPT__DT_USER;
   InputPara.Opt__DtLevel            = OPT__DT_LEVEL;
   InputPara.Opt__RecordDt           = OPT__RECORD_DT;
   InputPara.AutoReduceDt            = AUTO_REDUCE_DT;
   InputPara.AutoReduceDtFactor      = AUTO_REDUCE_DT_FACTOR;
   InputPara.AutoReduceDtFactorMin   = AUTO_REDUCE_DT_FACTOR_MIN;

// domain refinement
   InputPara.RegridCount             = REGRID_COUNT;
   InputPara.FlagBufferSize          = FLAG_BUFFER_SIZE;
   InputPara.FlagBufferSizeMaxM1Lv   = FLAG_BUFFER_SIZE_MAXM1_LV;
   InputPara.FlagBufferSizeMaxM2Lv   = FLAG_BUFFER_SIZE_MAXM2_LV;
   InputPara.MaxLevel                = MAX_LEVEL;
   InputPara.Opt__Flag_Rho           = OPT__FLAG_RHO;
   InputPara.Opt__Flag_RhoGradient   = OPT__FLAG_RHO_GRADIENT;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Flag_PresGradient  = OPT__FLAG_PRES_GRADIENT;
   InputPara.Opt__Flag_EngyGradient  = OPT__FLAG_ENGY_GRADIENT;
   InputPara.Opt__Flag_Vorticity     = OPT__FLAG_VORTICITY;
   InputPara.Opt__Flag_Jeans         = OPT__FLAG_JEANS;
#  ifdef MHD
   InputPara.Opt__Flag_Current       = OPT__FLAG_CURRENT;
#  endif
#  ifdef SRHD
   InputPara.Opt__Flag_4Velocity       = OPT__FLAG_4VELOCITY;
   InputPara.Opt__Flag_MomOverDens     = OPT__FLAG_MOM_OVER_DENS;
   InputPara.Opt__Flag_LorentzGradient = OPT__FLAG_LORENTZ_GRADIENT;
#  endif
#  endif
#  if ( MODEL == ELBDM )
   InputPara.Opt__Flag_EngyDensity   = OPT__FLAG_ENGY_DENSITY;
#  endif
   InputPara.Opt__Flag_LohnerDens    = OPT__FLAG_LOHNER_DENS;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Flag_LohnerEngy    = OPT__FLAG_LOHNER_ENGY;
   InputPara.Opt__Flag_LohnerPres    = OPT__FLAG_LOHNER_PRES;
   InputPara.Opt__Flag_LohnerTemp    = OPT__FLAG_LOHNER_TEMP;
#  endif
   InputPara.Opt__Flag_LohnerForm    = OPT__FLAG_LOHNER_FORM;
   InputPara.Opt__Flag_User          = OPT__FLAG_USER;
   InputPara.Opt__Flag_User_Num      = OPT__FLAG_USER_NUM;
   InputPara.Opt__Flag_Region        = OPT__FLAG_REGION;
#  ifdef PARTICLE
   InputPara.Opt__Flag_NParPatch     = OPT__FLAG_NPAR_PATCH;
   InputPara.Opt__Flag_NParCell      = OPT__FLAG_NPAR_CELL;
   InputPara.Opt__Flag_ParMassCell   = OPT__FLAG_PAR_MASS_CELL;
#  endif
   InputPara.Opt__NoFlagNearBoundary = OPT__NO_FLAG_NEAR_BOUNDARY;
   InputPara.Opt__PatchCount         = OPT__PATCH_COUNT;
#  ifdef PARTICLE
   InputPara.Opt__ParticleCount      = OPT__PARTICLE_COUNT;
#  endif
   InputPara.Opt__ReuseMemory        = OPT__REUSE_MEMORY;
   InputPara.Opt__MemoryPool         = OPT__MEMORY_POOL;

// load balance
#  ifdef LOAD_BALANCE
   InputPara.LB_WLI_Max              = amr->LB->WLI_Max;
#  ifdef PARTICLE
   InputPara.LB_Par_Weight           = amr->LB->Par_Weight;
#  endif
   InputPara.Opt__RecordLoadBalance  = OPT__RECORD_LOAD_BALANCE;
#  endif
   InputPara.Opt__MinimizeMPIBarrier = OPT__MINIMIZE_MPI_BARRIER;

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   InputPara.Gamma                   = GAMMA;
   InputPara.MolecularWeight         = MOLECULAR_WEIGHT;
   InputPara.IsoTemp                 = ISO_TEMP;
   InputPara.MinMod_Coeff            = MINMOD_COEFF;
   InputPara.Opt__LR_Limiter         = OPT__LR_LIMITER;
   InputPara.Opt__1stFluxCorr        = OPT__1ST_FLUX_CORR;
   InputPara.Opt__1stFluxCorrScheme  = OPT__1ST_FLUX_CORR_SCHEME;
#  endif

// ELBDM solvers
#  if ( MODEL == ELBDM )
   InputPara.ELBDM_Mass              = ELBDM_MASS;
   InputPara.ELBDM_PlanckConst       = ELBDM_PLANCK_CONST;
#  ifdef QUARTIC_SELF_INTERACTION
   InputPara.ELBDM_Lambda            = ELBDM_LAMBDA;
#  endif
   InputPara.ELBDM_Taylor3_Coeff     = ELBDM_TAYLOR3_COEFF;
   InputPara.ELBDM_Taylor3_Auto      = ELBDM_TAYLOR3_AUTO;
#  endif

// fluid solvers in different models
   InputPara.Flu_GPU_NPGroup         = FLU_GPU_NPGROUP;
   InputPara.GPU_NStream             = GPU_NSTREAM;
   InputPara.Opt__FixUp_Flux         = OPT__FIXUP_FLUX;
#  ifdef MHD
   InputPara.Opt__FixUp_Electric     = OPT__FIXUP_ELECTRIC;
#  endif
   InputPara.Opt__FixUp_Restrict     = OPT__FIXUP_RESTRICT;
   InputPara.Opt__CorrAfterAllSync   = OPT__CORR_AFTER_ALL_SYNC;
   InputPara.Opt__NormalizePassive   = OPT__NORMALIZE_PASSIVE;

   InputPara.NormalizePassive_NVar   = PassiveNorm_NVar;

   for (int v=0; v<NCOMP_PASSIVE; v++)
   InputPara.NormalizePassive_VarIdx[v] = PassiveNorm_VarIdx[v];

   for (int v=0; v<NCOMP_TOTAL; v++)
   InputPara.FieldLabel[v]           = FieldLabel[v];

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   InputPara.MagLabel[v]             = MagLabel[v];
#  endif

   InputPara.Opt__OverlapMPI         = OPT__OVERLAP_MPI;
   InputPara.Opt__ResetFluid         = OPT__RESET_FLUID;
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   InputPara.MinDens                 = MIN_DENS;
#  endif
#  if ( MODEL == HYDRO )
   InputPara.MinPres                 = MIN_PRES;
   InputPara.MinEint                 = MIN_EINT;
   InputPara.Opt__LastResortFloor    = OPT__LAST_RESORT_FLOOR;
   InputPara.JeansMinPres            = JEANS_MIN_PRES;
   InputPara.JeansMinPres_Level      = JEANS_MIN_PRES_LEVEL;
   InputPara.JeansMinPres_NCell      = JEANS_MIN_PRES_NCELL;
#  endif
#  ifdef DUAL_ENERGY
   InputPara.DualEnergySwitch        = DUAL_ENERGY_SWITCH;
#  endif

// self-gravity
#  ifdef GRAVITY
   InputPara.NewtonG                 = NEWTON_G;
#  if   ( POT_SCHEME == SOR )
   InputPara.SOR_Omega               = SOR_OMEGA;
   InputPara.SOR_MaxIter             = SOR_MAX_ITER;
   InputPara.SOR_MinIter             = SOR_MIN_ITER;
#  elif ( POT_SCHEME == MG )
   InputPara.MG_MaxIter              = MG_MAX_ITER;
   InputPara.MG_NPreSmooth           = MG_NPRE_SMOOTH;
   InputPara.MG_NPostSmooth          = MG_NPOST_SMOOTH;
   InputPara.MG_ToleratedError       = MG_TOLERATED_ERROR;
#  endif
   InputPara.Pot_GPU_NPGroup         = POT_GPU_NPGROUP;
   InputPara.Opt__GraP5Gradient      = OPT__GRA_P5_GRADIENT;
   InputPara.Opt__SelfGravity        = OPT__SELF_GRAVITY;
   InputPara.Opt__ExtAcc             = OPT__EXT_ACC;
   InputPara.Opt__ExtPot             = OPT__EXT_POT;
   InputPara.ExtPotTable_Name        = EXT_POT_TABLE_NAME;
   for (int d=0; d<3; d++)
   InputPara.ExtPotTable_NPoint[d]   = EXT_POT_TABLE_NPOINT[d];
   InputPara.ExtPotTable_dh          = EXT_POT_TABLE_DH;
   for (int d=0; d<3; d++)
   InputPara.ExtPotTable_EdgeL[d]    = EXT_POT_TABLE_EDGEL[d];
   InputPara.ExtPotTable_Float8      = EXT_POT_TABLE_FLOAT8;
   InputPara.Opt__GravityExtraMass   = OPT__GRAVITY_EXTRA_MASS;
#  endif

// source terms
   InputPara.Src_Deleptonization     = SrcTerms.Deleptonization;
   InputPara.Src_User                = SrcTerms.User;
   InputPara.Src_GPU_NPGroup         = SRC_GPU_NPGROUP;

// Grackle
#  ifdef SUPPORT_GRACKLE
   InputPara.Grackle_Activate        = GRACKLE_ACTIVATE;
   InputPara.Grackle_Verbose         = GRACKLE_VERBOSE;
   InputPara.Grackle_Cooling         = GRACKLE_COOLING;
   InputPara.Grackle_Primordial      = GRACKLE_PRIMORDIAL;
   InputPara.Grackle_Metal           = GRACKLE_METAL;
   InputPara.Grackle_UV              = GRACKLE_UV;
   InputPara.Grackle_CMB_Floor       = GRACKLE_CMB_FLOOR;
   InputPara.Grackle_PE_Heating      = GRACKLE_PE_HEATING;
   InputPara.Grackle_PE_HeatingRate  = GRACKLE_PE_HEATING_RATE;
   InputPara.Grackle_CloudyTable     = GRACKLE_CLOUDY_TABLE;
   InputPara.Grackle_ThreeBodyRate   = GRACKLE_THREE_BODY_RATE;
   InputPara.Grackle_CIE_Cooling     = GRACKLE_CIE_COOLING;
   InputPara.Grackle_H2_OpaApprox    = GRACKLE_H2_OPA_APPROX;
   InputPara.Che_GPU_NPGroup         = CHE_GPU_NPGROUP;
#  endif

// star formation
#  ifdef STAR_FORMATION
   InputPara.SF_CreateStar_Scheme       = SF_CREATE_STAR_SCHEME;
   InputPara.SF_CreateStar_RSeed        = SF_CREATE_STAR_RSEED;
   InputPara.SF_CreateStar_DetRandom    = SF_CREATE_STAR_DET_RANDOM;
   InputPara.SF_CreateStar_MinLevel     = SF_CREATE_STAR_MIN_LEVEL;
   InputPara.SF_CreateStar_MinGasDens   = SF_CREATE_STAR_MIN_GAS_DENS;
   InputPara.SF_CreateStar_MassEff      = SF_CREATE_STAR_MASS_EFF;
   InputPara.SF_CreateStar_MinStarMass  = SF_CREATE_STAR_MIN_STAR_MASS;
   InputPara.SF_CreateStar_MaxStarMFrac = SF_CREATE_STAR_MAX_STAR_MFRAC;
#  endif

// initialization
   InputPara.Opt__Init               = OPT__INIT;
   InputPara.RestartLoadNRank        = RESTART_LOAD_NRANK;
   InputPara.Opt__RestartReset       = OPT__RESTART_RESET;
   InputPara.Opt__UM_IC_Level        = OPT__UM_IC_LEVEL;
   InputPara.Opt__UM_IC_NVar         = OPT__UM_IC_NVAR;
   InputPara.Opt__UM_IC_Format       = OPT__UM_IC_FORMAT;
   InputPara.Opt__UM_IC_Downgrade    = OPT__UM_IC_DOWNGRADE;
   InputPara.Opt__UM_IC_Refine       = OPT__UM_IC_REFINE;
   InputPara.Opt__UM_IC_LoadNRank    = OPT__UM_IC_LOAD_NRANK;
   InputPara.Opt__InitRestrict       = OPT__INIT_RESTRICT;
   InputPara.Opt__InitGridWithOMP    = OPT__INIT_GRID_WITH_OMP;
   InputPara.Opt__GPUID_Select       = OPT__GPUID_SELECT;
   InputPara.Init_Subsampling_NCell  = INIT_SUBSAMPLING_NCELL;
#  ifdef MHD
   InputPara.Opt__InitBFieldByFile   = OPT__INIT_BFIELD_BYFILE;
#  endif

// interpolation schemes
   InputPara.Opt__Int_Time           = OPT__INT_TIME;
#  if ( MODEL == ELBDM )
   InputPara.Opt__Int_Phase          = OPT__INT_PHASE;
#  endif
   InputPara.Opt__Flu_IntScheme      = OPT__FLU_INT_SCHEME;
   InputPara.Opt__RefFlu_IntScheme   = OPT__REF_FLU_INT_SCHEME;
#  ifdef MHD
   InputPara.Opt__Mag_IntScheme      = OPT__MAG_INT_SCHEME;
   InputPara.Opt__RefMag_IntScheme   = OPT__REF_MAG_INT_SCHEME;
#  endif
#  ifdef GRAVITY
   InputPara.Opt__Pot_IntScheme      = OPT__POT_INT_SCHEME;
   InputPara.Opt__Rho_IntScheme      = OPT__RHO_INT_SCHEME;
   InputPara.Opt__Gra_IntScheme      = OPT__GRA_INT_SCHEME;
   InputPara.Opt__RefPot_IntScheme   = OPT__REF_POT_INT_SCHEME;
#  endif
   InputPara.IntMonoCoeff            = INT_MONO_COEFF;
   InputPara.IntOppSign0thOrder      = INT_OPP_SIGN_0TH_ORDER;

// data dump
   InputPara.Opt__Output_Total       = OPT__OUTPUT_TOTAL;
   InputPara.Opt__Output_Part        = OPT__OUTPUT_PART;
   InputPara.Opt__Output_User        = OPT__OUTPUT_USER;
#  ifdef PARTICLE
   InputPara.Opt__Output_ParText     = OPT__OUTPUT_PAR_TEXT;
#  endif
   InputPara.Opt__Output_BasePS      = OPT__OUTPUT_BASEPS;
   InputPara.Opt__Output_Base        = OPT__OUTPUT_BASE;
#  ifdef GRAVITY
   InputPara.Opt__Output_Pot         = OPT__OUTPUT_POT;
#  endif
#  ifdef PARTICLE
   InputPara.Opt__Output_ParDens     = OPT__OUTPUT_PAR_DENS;
#  endif
#  ifdef MHD
   InputPara.Opt__Output_CC_Mag      = OPT__OUTPUT_CC_MAG;
#  endif
   InputPara.Opt__Output_Mode        = OPT__OUTPUT_MODE;
   InputPara.Opt__Output_Step        = OUTPUT_STEP;
   InputPara.Opt__Output_Dt          = OUTPUT_DT;
   InputPara.Output_PartX            = OUTPUT_PART_X;
   InputPara.Output_PartY            = OUTPUT_PART_Y;
   InputPara.Output_PartZ            = OUTPUT_PART_Z;
   InputPara.InitDumpID              = INIT_DUMPID;

// miscellaneous
   InputPara.Opt__Verbose            = OPT__VERBOSE;
   InputPara.Opt__TimingBarrier      = OPT__TIMING_BARRIER;
   InputPara.Opt__TimingBalance      = OPT__TIMING_BALANCE;
   InputPara.Opt__TimingMPI          = OPT__TIMING_MPI;
   InputPara.Opt__RecordNote         = OPT__RECORD_NOTE;
   InputPara.Opt__RecordUnphy        = OPT__RECORD_UNPHY;
   InputPara.Opt__RecordMemory       = OPT__RECORD_MEMORY;
   InputPara.Opt__RecordPerformance  = OPT__RECORD_PERFORMANCE;
   InputPara.Opt__ManualControl      = OPT__MANUAL_CONTROL;
   InputPara.Opt__RecordUser         = OPT__RECORD_USER;
   InputPara.Opt__OptimizeAggressive = OPT__OPTIMIZE_AGGRESSIVE;

// simulation checks
   InputPara.Opt__Ck_Refine          = OPT__CK_REFINE;
   InputPara.Opt__Ck_ProperNesting   = OPT__CK_PROPER_NESTING;
   InputPara.Opt__Ck_Conservation    = OPT__CK_CONSERVATION;
   InputPara.Opt__Ck_NormPassive     = OPT__CK_NORMALIZE_PASSIVE;
   InputPara.Opt__Ck_Restrict        = OPT__CK_RESTRICT;
   InputPara.Opt__Ck_Finite          = OPT__CK_FINITE;
   InputPara.Opt__Ck_PatchAllocate   = OPT__CK_PATCH_ALLOCATE;
   InputPara.Opt__Ck_FluxAllocate    = OPT__CK_FLUX_ALLOCATE;
#  if ( MODEL == HYDRO )
   InputPara.Opt__Ck_Negative        = OPT__CK_NEGATIVE;
#  endif
   InputPara.Opt__Ck_MemFree         = OPT__CK_MEMFREE;
#  ifdef PARTICLE
   InputPara.Opt__Ck_Particle        = OPT__CK_PARTICLE;
#  endif
#  ifdef MHD
   InputPara.Opt__Ck_InterfaceB      = OPT__CK_INTERFACE_B;
   InputPara.Opt__Ck_DivergenceB     = OPT__CK_DIVERGENCE_B;
#  endif

// flag tables
#  if   ( MODEL == HYDRO )
   const bool Opt__FlagLohner = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP );
#  elif ( MODEL == ELBDM )
   const bool Opt__FlagLohner = OPT__FLAG_LOHNER_DENS;
#  endif

   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      InputPara.FlagTable_Rho         [lv]    = FlagTable_Rho         [lv];
      InputPara.FlagTable_RhoGradient [lv]    = FlagTable_RhoGradient [lv];

      for (int t=0; t<4; t++)
      InputPara.FlagTable_Lohner      [lv][t] = FlagTable_Lohner      [lv][t];

      InputPara.FlagTable_User        [lv].p   = malloc( OPT__FLAG_USER_NUM*sizeof(double) );
      InputPara.FlagTable_User        [lv].len = OPT__FLAG_USER_NUM;
      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      ( (double *) InputPara.FlagTable_User[lv].p )[t] = FlagTable_User[lv][t];

#     if   ( MODEL == HYDRO )
      InputPara.FlagTable_PresGradient[lv]    = FlagTable_PresGradient[lv];
      InputPara.FlagTable_EngyGradient[lv]    = FlagTable_EngyGradient[lv];
      InputPara.FlagTable_Vorticity   [lv]    = FlagTable_Vorticity   [lv];
      InputPara.FlagTable_Jeans       [lv]    = FlagTable_Jeans       [lv];
#     ifdef MHD
      InputPara.FlagTable_Current     [lv]    = FlagTable_Current     [lv];
#     endif
#     ifdef SRHD
      InputPara.FlagTable_LorentzFactorGradient[lv]    = FlagTable_LorentzFactorGradient[lv];
#     endif

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++)
      InputPara.FlagTable_EngyDensity [lv][t] = FlagTable_EngyDensity [lv][t];
#     endif

#     ifdef PARTICLE
      InputPara.FlagTable_NParPatch   [lv]    = FlagTable_NParPatch   [lv];
      InputPara.FlagTable_NParCell    [lv]    = FlagTable_NParCell    [lv];
      InputPara.FlagTable_ParMassCell [lv]    = FlagTable_ParMassCell [lv];
#     endif
   } // for (int lv=0; lv<NLEVEL-1; lv++)

} // FUNCTION : FillIn_InputPara



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_KeyInfo
// Description :  Create the HDF5 compound datatype for KeyInfo
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_KeyInfo( hid_t &H5_TypeID )
{

// create the array type
   const hsize_t H5_ArrDims_3Var         = 3;                        // array size of [3]
   const hsize_t H5_ArrDims_NLv          = NLEVEL;                   // array size of [NLEVEL]

   const hid_t   H5_TypeID_Arr_3Double   = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_3Var      );
   const hid_t   H5_TypeID_Arr_3Int      = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_3Var      );
   const hid_t   H5_TypeID_Arr_NLvInt    = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NLv       );
   const hid_t   H5_TypeID_Arr_NLvLong   = H5Tarray_create( H5T_NATIVE_LONG,   1, &H5_ArrDims_NLv       );
   const hid_t   H5_TypeID_Arr_NLvDouble = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_NLv       );


// create the "variable-length string" datatype
   hid_t  H5_TypeID_VarStr;
   herr_t H5_Status;

   H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status        = H5Tset_size( H5_TypeID_VarStr, H5T_VARIABLE );


// get the compound type
   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(KeyInfo_t) );

   H5Tinsert( H5_TypeID, "FormatVersion",        HOFFSET(KeyInfo_t,FormatVersion       ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Model",                HOFFSET(KeyInfo_t,Model               ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Float8",               HOFFSET(KeyInfo_t,Float8              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Gravity",              HOFFSET(KeyInfo_t,Gravity             ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "Particle",             HOFFSET(KeyInfo_t,Particle            ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NLevel",               HOFFSET(KeyInfo_t,NLevel              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NCompFluid",           HOFFSET(KeyInfo_t,NCompFluid          ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NCompPassive",         HOFFSET(KeyInfo_t,NCompPassive        ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "PatchSize",            HOFFSET(KeyInfo_t,PatchSize           ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "DumpID",               HOFFSET(KeyInfo_t,DumpID              ), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "NX0",                  HOFFSET(KeyInfo_t,NX0                 ), H5_TypeID_Arr_3Int      );
   H5Tinsert( H5_TypeID, "BoxScale",             HOFFSET(KeyInfo_t,BoxScale            ), H5_TypeID_Arr_3Int      );
   H5Tinsert( H5_TypeID, "NPatch",               HOFFSET(KeyInfo_t,NPatch              ), H5_TypeID_Arr_NLvInt    );
   H5Tinsert( H5_TypeID, "CellScale",            HOFFSET(KeyInfo_t,CellScale           ), H5_TypeID_Arr_NLvInt    );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Magnetohydrodynamics", HOFFSET(KeyInfo_t,Magnetohydrodynamics), H5T_NATIVE_INT          );
   H5Tinsert( H5_TypeID, "CosmicRay",            HOFFSET(KeyInfo_t,CosmicRay),            H5T_NATIVE_INT          );
#  endif

   H5Tinsert( H5_TypeID, "Step",                 HOFFSET(KeyInfo_t,Step                ), H5T_NATIVE_LONG         );
   H5Tinsert( H5_TypeID, "AdvanceCounter",       HOFFSET(KeyInfo_t,AdvanceCounter      ), H5_TypeID_Arr_NLvLong   );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_NPar",             HOFFSET(KeyInfo_t,Par_NPar),             H5T_NATIVE_LONG         );
   H5Tinsert( H5_TypeID, "Par_NAttStored",       HOFFSET(KeyInfo_t,Par_NAttStored      ), H5T_NATIVE_INT          );
#  endif

   H5Tinsert( H5_TypeID, "BoxSize",              HOFFSET(KeyInfo_t,BoxSize             ), H5_TypeID_Arr_3Double   );
   H5Tinsert( H5_TypeID, "Time",                 HOFFSET(KeyInfo_t,Time                ), H5_TypeID_Arr_NLvDouble );
   H5Tinsert( H5_TypeID, "CellSize",             HOFFSET(KeyInfo_t,CellSize            ), H5_TypeID_Arr_NLvDouble );
   H5Tinsert( H5_TypeID, "dTime_AllLv",          HOFFSET(KeyInfo_t,dTime_AllLv         ), H5_TypeID_Arr_NLvDouble );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "AveDens_Init",         HOFFSET(KeyInfo_t,AveDens_Init        ), H5T_NATIVE_DOUBLE       );
#  endif

   H5Tinsert( H5_TypeID, "CodeVersion",          HOFFSET(KeyInfo_t,CodeVersion         ), H5_TypeID_VarStr        );
   H5Tinsert( H5_TypeID, "DumpWallTime",         HOFFSET(KeyInfo_t,DumpWallTime        ), H5_TypeID_VarStr        );


// free memory
   H5_Status = H5Tclose( H5_TypeID_Arr_3Double       );
   H5_Status = H5Tclose( H5_TypeID_Arr_3Int          );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvInt        );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvLong       );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvDouble     );
   H5_Status = H5Tclose( H5_TypeID_VarStr            );

} // FUNCTION : GetCompound_KeyInfo



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_Makefile
// Description :  Create the HDF5 compound datatype for Makefile
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_Makefile( hid_t &H5_TypeID )
{

   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(Makefile_t) );

   H5Tinsert( H5_TypeID, "Model",                  HOFFSET(Makefile_t,Model                  ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Gravity",                HOFFSET(Makefile_t,Gravity                ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Comoving",               HOFFSET(Makefile_t,Comoving               ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Particle",               HOFFSET(Makefile_t,Particle               ), H5T_NATIVE_INT );

   H5Tinsert( H5_TypeID, "UseGPU",                 HOFFSET(Makefile_t,UseGPU                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "GAMER_Debug",            HOFFSET(Makefile_t,GAMER_Debug            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "BitwiseReproducibility", HOFFSET(Makefile_t,BitwiseReproducibility ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Timing",                 HOFFSET(Makefile_t,Timing                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "TimingSolver",           HOFFSET(Makefile_t,TimingSolver           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Float8",                 HOFFSET(Makefile_t,Float8                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Serial",                 HOFFSET(Makefile_t,Serial                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "LoadBalance",            HOFFSET(Makefile_t,LoadBalance            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "OverlapMPI",             HOFFSET(Makefile_t,OverlapMPI             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "OpenMP",                 HOFFSET(Makefile_t,OpenMP                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "GPU_Arch",               HOFFSET(Makefile_t,GPU_Arch               ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Laohu",                  HOFFSET(Makefile_t,Laohu                  ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportHDF5",            HOFFSET(Makefile_t,SupportHDF5            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportGSL",             HOFFSET(Makefile_t,SupportGSL             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SupportGrackle",         HOFFSET(Makefile_t,SupportGrackle         ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "RandomNumber",           HOFFSET(Makefile_t,RandomNumber           ), H5T_NATIVE_INT );

   H5Tinsert( H5_TypeID, "NLevel",                 HOFFSET(Makefile_t,NLevel                 ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "MaxPatch",               HOFFSET(Makefile_t,MaxPatch               ), H5T_NATIVE_INT );

#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "PotScheme",              HOFFSET(Makefile_t,PotScheme              ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "StorePotGhost",          HOFFSET(Makefile_t,StorePotGhost          ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "UnsplitGravity",         HOFFSET(Makefile_t,UnsplitGravity         ), H5T_NATIVE_INT );
#  endif

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "FluScheme",              HOFFSET(Makefile_t,FluScheme              ), H5T_NATIVE_INT );
#  ifdef LR_SCHEME
   H5Tinsert( H5_TypeID, "LRScheme",               HOFFSET(Makefile_t,LRScheme               ), H5T_NATIVE_INT );
#  endif
#  ifdef RSOLVER
   H5Tinsert( H5_TypeID, "RSolver",                HOFFSET(Makefile_t,RSolver                ), H5T_NATIVE_INT );
#  endif
   H5Tinsert( H5_TypeID, "DualEnergy",             HOFFSET(Makefile_t,DualEnergy             ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Magnetohydrodynamics",   HOFFSET(Makefile_t,Magnetohydrodynamics   ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SRHydrodynamics",        HOFFSET(Makefile_t,SRHydrodynamics        ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "CosmicRay",              HOFFSET(Makefile_t,CosmicRay              ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "EoS",                    HOFFSET(Makefile_t,EoS                    ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "BarotropicEoS",          HOFFSET(Makefile_t,BarotropicEoS          ), H5T_NATIVE_INT );

#  ifdef SRHD
   H5Tinsert( H5_TypeID, "SRHydroReducedEnergy",   HOFFSET(Makefile_t,SRHydroReducedEnergy  ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SRHydroFourVelocity",    HOFFSET(Makefile_t,SRHydroFourVelocity   ), H5T_NATIVE_INT );
#  endif

#  elif ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "ConserveMass",           HOFFSET(Makefile_t,ConserveMass           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Laplacian4th",           HOFFSET(Makefile_t,Laplacian4th           ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "SelfInteraction4",       HOFFSET(Makefile_t,SelfInteraction4       ), H5T_NATIVE_INT );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "StoreParAcc",            HOFFSET(Makefile_t,StoreParAcc            ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "StarFormation",          HOFFSET(Makefile_t,StarFormation          ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Par_NAttUser",           HOFFSET(Makefile_t,Par_NAttUser           ), H5T_NATIVE_INT );
#  endif

} // FUNCTION : GetCompound_Makefile



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_SymConst
// Description :  Create the HDF5 compound datatype for SymConst
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_SymConst( hid_t &H5_TypeID )
{

   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(SymConst_t) );

   H5Tinsert( H5_TypeID, "NCompFluid",           HOFFSET(SymConst_t,NCompFluid          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NCompPassive",         HOFFSET(SymConst_t,NCompPassive        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "PatchSize",            HOFFSET(SymConst_t,PatchSize           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn",              HOFFSET(SymConst_t,Flu_NIn             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NOut",             HOFFSET(SymConst_t,Flu_NOut            ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn_T",            HOFFSET(SymConst_t,Flu_NIn_T           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NIn_S",            HOFFSET(SymConst_t,Flu_NIn_S           ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_NOut_S",           HOFFSET(SymConst_t,Flu_NOut_S          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NFluxFluid",           HOFFSET(SymConst_t,NFluxFluid          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "NFluxPassive",         HOFFSET(SymConst_t,NFluxPassive        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_GhostSize",        HOFFSET(SymConst_t,Flu_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_Nxt",              HOFFSET(SymConst_t,Flu_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Debug_HDF5",           HOFFSET(SymConst_t,Debug_HDF5          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SibOffsetNonperiodic", HOFFSET(SymConst_t,SibOffsetNonperiodic), H5T_NATIVE_INT    );
#  ifdef LOAD_BALANCE
   H5Tinsert( H5_TypeID, "SonOffsetLB",          HOFFSET(SymConst_t,SonOffsetLB         ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "TinyNumber",           HOFFSET(SymConst_t,TinyNumber          ), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "HugeNumber",           HOFFSET(SymConst_t,HugeNumber          ), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "MaxError",             HOFFSET(SymConst_t,MaxError            ), H5T_NATIVE_DOUBLE );

#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Gra_NIn",              HOFFSET(SymConst_t,Gra_NIn             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Pot_GhostSize",        HOFFSET(SymConst_t,Pot_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_GhostSize",        HOFFSET(SymConst_t,Gra_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Rho_GhostSize",        HOFFSET(SymConst_t,Rho_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Pot_Nxt",              HOFFSET(SymConst_t,Pot_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_Nxt",              HOFFSET(SymConst_t,Gra_Nxt             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Rho_Nxt",              HOFFSET(SymConst_t,Rho_Nxt             ), H5T_NATIVE_INT    );
#  ifdef UNSPLIT_GRAVITY
   H5Tinsert( H5_TypeID, "USG_GhostSizeF",       HOFFSET(SymConst_t,USG_GhostSizeF      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_GhostSizeG",       HOFFSET(SymConst_t,USG_GhostSizeG      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_NxtF",             HOFFSET(SymConst_t,USG_NxtF            ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "USG_NxtG",             HOFFSET(SymConst_t,USG_NxtG            ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "ExtPot_BlockSize",     HOFFSET(SymConst_t,ExtPot_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Gra_BlockSize",        HOFFSET(SymConst_t,Gra_BlockSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ExtPotNAuxMax",        HOFFSET(SymConst_t,ExtPotNAuxMax       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ExtAccNAuxMax",        HOFFSET(SymConst_t,ExtAccNAuxMax       ), H5T_NATIVE_INT    );
#  if   ( POT_SCHEME == SOR )
   H5Tinsert( H5_TypeID, "Pot_BlockSize_z",      HOFFSET(SymConst_t,Pot_BlockSize_z     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "UsePSolver_10to14",    HOFFSET(SymConst_t,UsePSolver_10to14   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_RhoShared",        HOFFSET(SymConst_t,SOR_RhoShared       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_CPotShared",       HOFFSET(SymConst_t,SOR_CPotShared      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_UseShuffle",       HOFFSET(SymConst_t,SOR_UseShuffle      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_UsePadding",       HOFFSET(SymConst_t,SOR_UsePadding      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "SOR_ModReduction",     HOFFSET(SymConst_t,SOR_ModReduction    ), H5T_NATIVE_INT    );
#  elif ( POT_SCHEME == MG  )
   H5Tinsert( H5_TypeID, "Pot_BlockSize_x",      HOFFSET(SymConst_t,Pot_BlockSize_x     ), H5T_NATIVE_INT    );
#  endif
#  endif // #ifdef GRAVITY

#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_NAttStored",       HOFFSET(SymConst_t,Par_NAttStored      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "RhoExt_GhostSize",     HOFFSET(SymConst_t,RhoExt_GhostSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Debug_Particle",       HOFFSET(SymConst_t,Debug_Particle      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "ParList_GrowthFactor", HOFFSET(SymConst_t,ParList_GrowthFactor), H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_TypeID, "ParList_ReduceFactor", HOFFSET(SymConst_t,ParList_ReduceFactor), H5T_NATIVE_DOUBLE );
#  endif

   H5Tinsert( H5_TypeID, "BitRep_Flux",          HOFFSET(SymConst_t,BitRep_Flux         ), H5T_NATIVE_INT    );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "BitRep_Electric",      HOFFSET(SymConst_t,BitRep_Electric     ), H5T_NATIVE_INT    );
#  endif

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Flu_BlockSize_x",      HOFFSET(SymConst_t,Flu_BlockSize_x     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_BlockSize_y",      HOFFSET(SymConst_t,Flu_BlockSize_y     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CheckNegativeInFluid", HOFFSET(SymConst_t,CheckNegativeInFluid), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CharReconstruction",   HOFFSET(SymConst_t,CharReconstruction  ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "LR_Eint",              HOFFSET(SymConst_t,LR_Eint             ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "CheckIntermediate",    HOFFSET(SymConst_t,CheckIntermediate   ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLL_NoRefState",       HOFFSET(SymConst_t,HLL_NoRefState      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLL_IncludeAllWaves",  HOFFSET(SymConst_t,HLL_IncludeAllWaves ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLLC_WaveSpeed",       HOFFSET(SymConst_t,HLLC_WaveSpeed      ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "HLLE_WaveSpeed",       HOFFSET(SymConst_t,HLLE_WaveSpeed      ), H5T_NATIVE_INT    );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "HLLD_WaveSpeed",       HOFFSET(SymConst_t,HLLD_WaveSpeed      ), H5T_NATIVE_INT    );
#  endif
#  ifdef N_FC_VAR
   H5Tinsert( H5_TypeID, "N_FC_Var",             HOFFSET(SymConst_t,N_FC_Var            ), H5T_NATIVE_INT    );
#  endif
#  ifdef N_SLOPE_PPM
   H5Tinsert( H5_TypeID, "N_Slope_PPM",          HOFFSET(SymConst_t,N_Slope_PPM         ), H5T_NATIVE_INT    );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "EulerY",               HOFFSET(SymConst_t,EulerY              ), H5T_NATIVE_INT    );
#  endif
   H5Tinsert( H5_TypeID, "EoSNAuxMax",           HOFFSET(SymConst_t,EoSNAuxMax          ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "EoSNTableMax",         HOFFSET(SymConst_t,EoSNTableMax        ), H5T_NATIVE_INT    );

#  elif  ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Flu_BlockSize_x",      HOFFSET(SymConst_t,Flu_BlockSize_x     ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Flu_BlockSize_y",      HOFFSET(SymConst_t,Flu_BlockSize_y     ), H5T_NATIVE_INT    );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   H5Tinsert( H5_TypeID, "dt_Flu_BlockSize",     HOFFSET(SymConst_t,dt_Flu_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "dt_Flu_UseShuffle",    HOFFSET(SymConst_t,dt_Flu_UseShuffle   ), H5T_NATIVE_INT    );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "dt_Gra_BlockSize",     HOFFSET(SymConst_t,dt_Gra_BlockSize    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "dt_Gra_UseShuffle",    HOFFSET(SymConst_t,dt_Gra_UseShuffle   ), H5T_NATIVE_INT    );
#  endif

   H5Tinsert( H5_TypeID, "Src_BlockSize",        HOFFSET(SymConst_t,Src_BlockSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_GhostSize",        HOFFSET(SymConst_t,Src_GhostSize       ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_NAuxDlep",         HOFFSET(SymConst_t,Src_NAuxDlep        ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_DlepProfNVar",     HOFFSET(SymConst_t,Src_DlepProfNVar    ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_DlepProfNBinMax",  HOFFSET(SymConst_t,Src_DlepProfNBinMax ), H5T_NATIVE_INT    );
   H5Tinsert( H5_TypeID, "Src_NAuxUser",         HOFFSET(SymConst_t,Src_NAuxUser        ), H5T_NATIVE_INT    );

} // FUNCTION : GetCompound_SymConst



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCompound_InputPara
// Description :  Create the HDF5 compound datatype for InputPara
//
// Note        :  1. Data structure is defined in "HDF5_Typedef.h"
//                2. The returned H5_TypeID must be closed manually
//                3. Call-by-reference
//
// Parameter   :  H5_TypeID : HDF5 type ID for storing the compound datatype
//-------------------------------------------------------------------------------------------------------
void GetCompound_InputPara( hid_t &H5_TypeID )
{

// create the array type
   const hsize_t H5_ArrDims_3Var              = 3;                    // array size of [3]
   const hsize_t H5_ArrDims_6Var              = 6;                    // array size of [6]
#  if ( NCOMP_PASSIVE > 0 )
   const hsize_t H5_ArrDims_NPassive          = NCOMP_PASSIVE;        // array size of [NCOMP_PASSIVE]
#  endif
#  if ( NLEVEL > 1 )
   const hsize_t H5_ArrDims_NLvM1             = NLEVEL-1;             // array size of [NLEVEL-1]
   const hsize_t H5_ArrDims_NLvM1_2[2]        = { NLEVEL-1, 2 };      // array size of [NLEVEL-1][2]
   const hsize_t H5_ArrDims_NLvM1_4[2]        = { NLEVEL-1, 4 };      // array size of [NLEVEL-1][4]
#  endif

   const hid_t   H5_TypeID_Arr_3Int           = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_3Var      );
   const hid_t   H5_TypeID_Arr_6Int           = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_6Var      );
#  if ( NCOMP_PASSIVE > 0 )
   const hid_t   H5_TypeID_Arr_NPassive       = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NPassive  );
#  endif
#  if ( NLEVEL > 1 )
   const hid_t   H5_TypeID_Arr_NLvM1Int       = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_NLvM1     );
   const hid_t   H5_TypeID_Arr_NLvM1Double    = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_NLvM1     );
   const hid_t   H5_TypeID_Arr_NLvM1_2Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_2   );
   const hid_t   H5_TypeID_Arr_NLvM1_4Double  = H5Tarray_create( H5T_NATIVE_DOUBLE, 2,  H5_ArrDims_NLvM1_4   );
   const hid_t   H5_TypeID_Arr_NLvM1_VLDouble = H5Tvlen_create ( H5T_NATIVE_DOUBLE );
#  endif
   const hid_t   H5_TypeID_Arr_3Double        = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_3Var      );


// create the "variable-length string" datatype
   hid_t  H5_TypeID_VarStr;
   herr_t H5_Status;

   H5_TypeID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status        = H5Tset_size( H5_TypeID_VarStr, H5T_VARIABLE );


// get the size of a single pointer, which is used for storing the array of variable-length strings
// --> FieldLabel[], MagLabel[], ParAttLabel[]
   const int PtrSize     = sizeof( char* );
   const int PtrSize_hvl = sizeof( hvl_t );
   char Key[MAX_STRING];


// get the compound type
   H5_TypeID = H5Tcreate( H5T_COMPOUND, sizeof(InputPara_t) );

// simulation scale
   H5Tinsert( H5_TypeID, "BoxSize",                 HOFFSET(InputPara_t,BoxSize                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "NX0_Tot",                 HOFFSET(InputPara_t,NX0_Tot                ), H5_TypeID_Arr_3Int );
   H5Tinsert( H5_TypeID, "MPI_NRank",               HOFFSET(InputPara_t,MPI_NRank              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "MPI_NRank_X",             HOFFSET(InputPara_t,MPI_NRank_X            ), H5_TypeID_Arr_3Int );
   H5Tinsert( H5_TypeID, "OMP_NThread",             HOFFSET(InputPara_t,OMP_NThread            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "EndT",                    HOFFSET(InputPara_t,EndT                   ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "EndStep",                 HOFFSET(InputPara_t,EndStep                ), H5T_NATIVE_LONG    );

// test problems
   H5Tinsert( H5_TypeID, "TestProb_ID",             HOFFSET(InputPara_t,TestProb_ID            ), H5T_NATIVE_INT     );

// code units
   H5Tinsert( H5_TypeID, "Opt__Unit",               HOFFSET(InputPara_t,Opt__Unit              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Unit_L",                  HOFFSET(InputPara_t,Unit_L                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_M",                  HOFFSET(InputPara_t,Unit_M                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_T",                  HOFFSET(InputPara_t,Unit_T                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_V",                  HOFFSET(InputPara_t,Unit_V                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_D",                  HOFFSET(InputPara_t,Unit_D                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_E",                  HOFFSET(InputPara_t,Unit_E                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Unit_P",                  HOFFSET(InputPara_t,Unit_P                 ), H5T_NATIVE_DOUBLE  );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Unit_B",                  HOFFSET(InputPara_t,Unit_B                 ), H5T_NATIVE_DOUBLE  );
#  endif

// boundary condition
   H5Tinsert( H5_TypeID, "Opt__BC_Flu",             HOFFSET(InputPara_t,Opt__BC_Flu            ), H5_TypeID_Arr_6Int );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__BC_Pot",             HOFFSET(InputPara_t,Opt__BC_Pot            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "GFunc_Coeff0",            HOFFSET(InputPara_t,GFunc_Coeff0           ), H5T_NATIVE_DOUBLE  );
#  endif

// particle
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Par_Init",                HOFFSET(InputPara_t,Par_Init               ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICFormat",            HOFFSET(InputPara_t,Par_ICFormat           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ICMass",              HOFFSET(InputPara_t,Par_ICMass             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Par_Interp",              HOFFSET(InputPara_t,Par_Interp             ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_Integ",               HOFFSET(InputPara_t,Par_Integ              ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_ImproveAcc",          HOFFSET(InputPara_t,Par_ImproveAcc         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_PredictPos",          HOFFSET(InputPara_t,Par_PredictPos         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Par_RemoveCell",          HOFFSET(InputPara_t,Par_RemoveCell         ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Par_GhostSize",           HOFFSET(InputPara_t,Par_GhostSize          ), H5T_NATIVE_INT     );

// store the name of all particle attributes
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
//    key for each particle attribute
      sprintf( Key, "ParAttLabel%02d", v );

//    assuming the offset between successive ParAttLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,ParAttLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }
#  endif

// cosmology
#  ifdef COMOVING
   H5Tinsert( H5_TypeID, "A_Init",                  HOFFSET(InputPara_t,A_Init                 ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "OmegaM0",                 HOFFSET(InputPara_t,OmegaM0                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Hubble0",                 HOFFSET(InputPara_t,Hubble0                ), H5T_NATIVE_DOUBLE  );
#  endif

// time-step determination
   H5Tinsert( H5_TypeID, "Dt__Max",                 HOFFSET(InputPara_t,Dt__Max                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__Fluid",               HOFFSET(InputPara_t,Dt__Fluid              ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__FluidInit",           HOFFSET(InputPara_t,Dt__FluidInit          ), H5T_NATIVE_DOUBLE  );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Dt__Gravity",             HOFFSET(InputPara_t,Dt__Gravity            ), H5T_NATIVE_DOUBLE  );
#  endif
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Dt__Phase",               HOFFSET(InputPara_t,Dt__Phase              ), H5T_NATIVE_DOUBLE  );
#  endif
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Dt__ParVel",              HOFFSET(InputPara_t,Dt__ParVel             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__ParVelMax",           HOFFSET(InputPara_t,Dt__ParVelMax          ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__ParAcc",              HOFFSET(InputPara_t,Dt__ParAcc             ), H5T_NATIVE_DOUBLE  );
#  endif
#  ifdef COMOVING
   H5Tinsert( H5_TypeID, "Dt__MaxDeltaA",           HOFFSET(InputPara_t,Dt__MaxDeltaA          ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "Dt__SyncParentLv",        HOFFSET(InputPara_t,Dt__SyncParentLv       ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Dt__SyncChildrenLv",      HOFFSET(InputPara_t,Dt__SyncChildrenLv     ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Opt__DtUser",             HOFFSET(InputPara_t,Opt__DtUser            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__DtLevel",            HOFFSET(InputPara_t,Opt__DtLevel           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__RecordDt",           HOFFSET(InputPara_t,Opt__RecordDt          ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "AutoReduceDt",            HOFFSET(InputPara_t,AutoReduceDt           ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "AutoReduceDtFactor",      HOFFSET(InputPara_t,AutoReduceDtFactor     ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "AutoReduceDtFactorMin",   HOFFSET(InputPara_t,AutoReduceDtFactorMin  ), H5T_NATIVE_DOUBLE  );


// domain refinement
   H5Tinsert( H5_TypeID, "RegridCount",             HOFFSET(InputPara_t,RegridCount            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSize",          HOFFSET(InputPara_t,FlagBufferSize         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSizeMaxM1Lv",   HOFFSET(InputPara_t,FlagBufferSizeMaxM1Lv  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "FlagBufferSizeMaxM2Lv",   HOFFSET(InputPara_t,FlagBufferSizeMaxM2Lv  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "MaxLevel",                HOFFSET(InputPara_t,MaxLevel               ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Rho",           HOFFSET(InputPara_t,Opt__Flag_Rho          ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_RhoGradient",   HOFFSET(InputPara_t,Opt__Flag_RhoGradient  ), H5T_NATIVE_INT     );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Flag_PresGradient",  HOFFSET(InputPara_t,Opt__Flag_PresGradient ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_EngyGradient",  HOFFSET(InputPara_t,Opt__Flag_EngyGradient ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Vorticity",     HOFFSET(InputPara_t,Opt__Flag_Vorticity    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Jeans",         HOFFSET(InputPara_t,Opt__Flag_Jeans        ), H5T_NATIVE_INT     );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Flag_Current",       HOFFSET(InputPara_t,Opt__Flag_Current      ), H5T_NATIVE_INT     );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "Opt__Flag_4Velocity",     HOFFSET(InputPara_t,Opt__Flag_4Velocity      ), H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Opt__Flag_MomOverDens",   HOFFSET(InputPara_t,Opt__Flag_MomOverDens  ),   H5T_NATIVE_INT );
   H5Tinsert( H5_TypeID, "Opt__Flag_LorentzGradient", HOFFSET(InputPara_t,Opt__Flag_LorentzGradient), H5T_NATIVE_INT );
#  endif
#  endif
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Opt__Flag_EngyDensity",   HOFFSET(InputPara_t,Opt__Flag_EngyDensity  ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerDens",    HOFFSET(InputPara_t,Opt__Flag_LohnerDens   ), H5T_NATIVE_INT     );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerEngy",    HOFFSET(InputPara_t,Opt__Flag_LohnerEngy   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerPres",    HOFFSET(InputPara_t,Opt__Flag_LohnerPres   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerTemp",    HOFFSET(InputPara_t,Opt__Flag_LohnerTemp   ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Flag_LohnerForm",    HOFFSET(InputPara_t,Opt__Flag_LohnerForm   ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_User",          HOFFSET(InputPara_t,Opt__Flag_User         ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_User_Num",      HOFFSET(InputPara_t,Opt__Flag_User_Num     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_Region",        HOFFSET(InputPara_t,Opt__Flag_Region       ), H5T_NATIVE_INT     );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Flag_NParPatch",     HOFFSET(InputPara_t,Opt__Flag_NParPatch    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_NParCell",      HOFFSET(InputPara_t,Opt__Flag_NParCell     ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__Flag_ParMassCell",   HOFFSET(InputPara_t,Opt__Flag_ParMassCell  ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__NoFlagNearBoundary", HOFFSET(InputPara_t,Opt__NoFlagNearBoundary), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__PatchCount",         HOFFSET(InputPara_t,Opt__PatchCount        ), H5T_NATIVE_INT     );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__ParticleCount",      HOFFSET(InputPara_t,Opt__ParticleCount     ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__ReuseMemory",        HOFFSET(InputPara_t,Opt__ReuseMemory       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__MemoryPool",         HOFFSET(InputPara_t,Opt__MemoryPool        ), H5T_NATIVE_INT     );

// load balance
#  ifdef LOAD_BALANCE
   H5Tinsert( H5_TypeID, "LB_WLI_Max",              HOFFSET(InputPara_t,LB_WLI_Max             ), H5T_NATIVE_DOUBLE  );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "LB_Par_Weight",           HOFFSET(InputPara_t,LB_Par_Weight          ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "Opt__RecordLoadBalance",  HOFFSET(InputPara_t,Opt__RecordLoadBalance ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__MinimizeMPIBarrier", HOFFSET(InputPara_t,Opt__MinimizeMPIBarrier), H5T_NATIVE_INT     );

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Gamma",                   HOFFSET(InputPara_t,Gamma                  ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MolecularWeight",         HOFFSET(InputPara_t,MolecularWeight        ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "IsoTemp",                 HOFFSET(InputPara_t,IsoTemp                ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "MinMod_Coeff",            HOFFSET(InputPara_t,MinMod_Coeff           ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "Opt__LR_Limiter",         HOFFSET(InputPara_t,Opt__LR_Limiter        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__1stFluxCorr",        HOFFSET(InputPara_t,Opt__1stFluxCorr       ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__1stFluxCorrScheme",  HOFFSET(InputPara_t,Opt__1stFluxCorrScheme ), H5T_NATIVE_INT     );
#  endif

// ELBDM solvers
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "ELBDM_Mass",              HOFFSET(InputPara_t,ELBDM_Mass             ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "ELBDM_PlanckConst",       HOFFSET(InputPara_t,ELBDM_PlanckConst      ), H5T_NATIVE_DOUBLE  );
#  ifdef QUARTIC_SELF_INTERACTION
   H5Tinsert( H5_TypeID, "ELBDM_Lambda",            HOFFSET(InputPara_t,ELBDM_Lambda           ), H5T_NATIVE_DOUBLE  );
#  endif
   H5Tinsert( H5_TypeID, "ELBDM_Taylor3_Coeff",     HOFFSET(InputPara_t,ELBDM_Taylor3_Coeff    ), H5T_NATIVE_DOUBLE  );
   H5Tinsert( H5_TypeID, "ELBDM_Taylor3_Auto",      HOFFSET(InputPara_t,ELBDM_Taylor3_Auto     ), H5T_NATIVE_INT     );
#  endif

// fluid solvers in different models
   H5Tinsert( H5_TypeID, "Flu_GPU_NPGroup",         HOFFSET(InputPara_t,Flu_GPU_NPGroup        ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "GPU_NStream",             HOFFSET(InputPara_t,GPU_NStream            ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__FixUp_Flux",         HOFFSET(InputPara_t,Opt__FixUp_Flux        ), H5T_NATIVE_INT     );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__FixUp_Electric",     HOFFSET(InputPara_t,Opt__FixUp_Electric    ), H5T_NATIVE_INT     );
#  endif
   H5Tinsert( H5_TypeID, "Opt__FixUp_Restrict",     HOFFSET(InputPara_t,Opt__FixUp_Restrict    ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__CorrAfterAllSync",   HOFFSET(InputPara_t,Opt__CorrAfterAllSync  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "Opt__NormalizePassive",   HOFFSET(InputPara_t,Opt__NormalizePassive  ), H5T_NATIVE_INT     );
   H5Tinsert( H5_TypeID, "NormalizePassive_NVar",   HOFFSET(InputPara_t,NormalizePassive_NVar  ), H5T_NATIVE_INT     );
#  if ( NCOMP_PASSIVE > 0 )
   H5Tinsert( H5_TypeID, "NormalizePassive_VarIdx", HOFFSET(InputPara_t,NormalizePassive_VarIdx), H5_TypeID_Arr_NPassive );
#  endif

// store the name of all fields
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
//    key for each field
      sprintf( Key, "FieldLabel%02d", v );

//    assuming the offset between successive FieldLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,FieldLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }

#  ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++)
   {
//    key for each field
      sprintf( Key, "MagLabel%02d", v );

//    assuming the offset between successive MagLabel pointers is "PtrSize", which is equal to "sizeof( char* )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,MagLabel[0])+v*PtrSize, H5_TypeID_VarStr );
   }
#  endif

   H5Tinsert( H5_TypeID, "Opt__OverlapMPI",         HOFFSET(InputPara_t,Opt__OverlapMPI        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ResetFluid",         HOFFSET(InputPara_t,Opt__ResetFluid        ), H5T_NATIVE_INT              );
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "MinDens",                 HOFFSET(InputPara_t,MinDens                ), H5T_NATIVE_DOUBLE           );
#  endif
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "MinPres",                 HOFFSET(InputPara_t,MinPres                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "MinEint",                 HOFFSET(InputPara_t,MinEint                ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Opt__LastResortFloor",    HOFFSET(InputPara_t,Opt__LastResortFloor   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres",            HOFFSET(InputPara_t,JeansMinPres           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres_Level",      HOFFSET(InputPara_t,JeansMinPres_Level     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "JeansMinPres_NCell",      HOFFSET(InputPara_t,JeansMinPres_NCell     ), H5T_NATIVE_INT              );
#  endif
#  ifdef DUAL_ENERGY
   H5Tinsert( H5_TypeID, "DualEnergySwitch",        HOFFSET(InputPara_t,DualEnergySwitch       ), H5T_NATIVE_DOUBLE           );
#  endif

// self-gravity
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "NewtonG",                 HOFFSET(InputPara_t,NewtonG                ), H5T_NATIVE_DOUBLE           );
#  if   ( POT_SCHEME == SOR )
   H5Tinsert( H5_TypeID, "SOR_Omega",               HOFFSET(InputPara_t,SOR_Omega              ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "SOR_MaxIter",             HOFFSET(InputPara_t,SOR_MaxIter            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "SOR_MinIter",             HOFFSET(InputPara_t,SOR_MinIter            ), H5T_NATIVE_INT              );
#  elif ( POT_SCHEME == MG )
   H5Tinsert( H5_TypeID, "MG_MaxIter",              HOFFSET(InputPara_t,MG_MaxIter             ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_NPreSmooth",           HOFFSET(InputPara_t,MG_NPreSmooth          ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_NPostSmooth",          HOFFSET(InputPara_t,MG_NPostSmooth         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "MG_ToleratedError",       HOFFSET(InputPara_t,MG_ToleratedError      ), H5T_NATIVE_DOUBLE           );
#  endif
   H5Tinsert( H5_TypeID, "Pot_GPU_NPGroup",         HOFFSET(InputPara_t,Pot_GPU_NPGroup        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GraP5Gradient",      HOFFSET(InputPara_t,Opt__GraP5Gradient     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__SelfGravity",        HOFFSET(InputPara_t,Opt__SelfGravity       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ExtAcc",             HOFFSET(InputPara_t,Opt__ExtAcc            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ExtPot",             HOFFSET(InputPara_t,Opt__ExtPot            ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "ExtPotTable_Name",        HOFFSET(InputPara_t,ExtPotTable_Name       ), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "ExtPotTable_NPoint",      HOFFSET(InputPara_t,ExtPotTable_NPoint     ), H5_TypeID_Arr_3Int          );
   H5Tinsert( H5_TypeID, "ExtPotTable_dh",          HOFFSET(InputPara_t,ExtPotTable_dh         ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "ExtPotTable_EdgeL",       HOFFSET(InputPara_t,ExtPotTable_EdgeL      ), H5_TypeID_Arr_3Double       );
   H5Tinsert( H5_TypeID, "ExtPotTable_Float8",      HOFFSET(InputPara_t,ExtPotTable_Float8     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GravityExtraMass",   HOFFSET(InputPara_t,Opt__GravityExtraMass  ), H5T_NATIVE_INT              );
#  endif // #ifdef GRAVITY

// source terms
   H5Tinsert( H5_TypeID, "Src_Deleptonization",     HOFFSET(InputPara_t,Src_Deleptonization    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Src_User",                HOFFSET(InputPara_t,Src_User               ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Src_GPU_NPGroup",         HOFFSET(InputPara_t,Src_GPU_NPGroup        ), H5T_NATIVE_INT              );

// Grackle
#  ifdef SUPPORT_GRACKLE
   H5Tinsert( H5_TypeID, "Grackle_Activate",        HOFFSET(InputPara_t,Grackle_Activate       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Verbose",         HOFFSET(InputPara_t,Grackle_Verbose        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Cooling",         HOFFSET(InputPara_t,Grackle_Cooling        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Primordial",      HOFFSET(InputPara_t,Grackle_Primordial     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_Metal",           HOFFSET(InputPara_t,Grackle_Metal          ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_UV",              HOFFSET(InputPara_t,Grackle_UV             ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_CMB_Floor",       HOFFSET(InputPara_t,Grackle_CMB_Floor      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_PE_Heating",      HOFFSET(InputPara_t,Grackle_PE_Heating     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_PE_HeatingRate",  HOFFSET(InputPara_t,Grackle_PE_HeatingRate ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Grackle_CloudyTable",     HOFFSET(InputPara_t,Grackle_CloudyTable    ), H5_TypeID_VarStr            );
   H5Tinsert( H5_TypeID, "Grackle_ThreeBodyRate",   HOFFSET(InputPara_t,Grackle_ThreeBodyRate  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_CIE_Cooling",     HOFFSET(InputPara_t,Grackle_CIE_Cooling    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Grackle_H2_OpaApprox",    HOFFSET(InputPara_t,Grackle_H2_OpaApprox   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Che_GPU_NPGroup",         HOFFSET(InputPara_t,Che_GPU_NPGroup        ), H5T_NATIVE_INT              );
#  endif

// star formation
#  ifdef STAR_FORMATION
   H5Tinsert( H5_TypeID, "SF_CreateStar_Scheme",       HOFFSET(InputPara_t,SF_CreateStar_Scheme       ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_RSeed",        HOFFSET(InputPara_t,SF_CreateStar_RSeed        ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_DetRandom",    HOFFSET(InputPara_t,SF_CreateStar_DetRandom    ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinLevel",     HOFFSET(InputPara_t,SF_CreateStar_MinLevel     ), H5T_NATIVE_INT       );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinGasDens",   HOFFSET(InputPara_t,SF_CreateStar_MinGasDens   ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MassEff",      HOFFSET(InputPara_t,SF_CreateStar_MassEff      ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MinStarMass",  HOFFSET(InputPara_t,SF_CreateStar_MinStarMass  ), H5T_NATIVE_DOUBLE    );
   H5Tinsert( H5_TypeID, "SF_CreateStar_MaxStarMFrac", HOFFSET(InputPara_t,SF_CreateStar_MaxStarMFrac ), H5T_NATIVE_DOUBLE    );
#  endif

// initialization
   H5Tinsert( H5_TypeID, "Opt__Init",               HOFFSET(InputPara_t,Opt__Init              ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "RestartLoadNRank",        HOFFSET(InputPara_t,RestartLoadNRank       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RestartReset",       HOFFSET(InputPara_t,Opt__RestartReset      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Level",        HOFFSET(InputPara_t,Opt__UM_IC_Level       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_NVar",         HOFFSET(InputPara_t,Opt__UM_IC_NVar        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Format",       HOFFSET(InputPara_t,Opt__UM_IC_Format      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Downgrade",    HOFFSET(InputPara_t,Opt__UM_IC_Downgrade   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_Refine",       HOFFSET(InputPara_t,Opt__UM_IC_Refine      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__UM_IC_LoadNRank",    HOFFSET(InputPara_t,Opt__UM_IC_LoadNRank   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__InitRestrict",       HOFFSET(InputPara_t,Opt__InitRestrict      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__InitGridWithOMP",    HOFFSET(InputPara_t,Opt__InitGridWithOMP   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__GPUID_Select",       HOFFSET(InputPara_t,Opt__GPUID_Select      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Init_Subsampling_NCell",  HOFFSET(InputPara_t,Init_Subsampling_NCell ), H5T_NATIVE_INT              );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__InitBFieldByFile",   HOFFSET(InputPara_t,Opt__InitBFieldByFile  ), H5T_NATIVE_INT              );
#  endif

// interpolation schemes
   H5Tinsert( H5_TypeID, "Opt__Int_Time",           HOFFSET(InputPara_t,Opt__Int_Time          ), H5T_NATIVE_INT              );
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "Opt__Int_Phase",          HOFFSET(InputPara_t,Opt__Int_Phase         ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Flu_IntScheme",      HOFFSET(InputPara_t,Opt__Flu_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefFlu_IntScheme",   HOFFSET(InputPara_t,Opt__RefFlu_IntScheme  ), H5T_NATIVE_INT              );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Mag_IntScheme",      HOFFSET(InputPara_t,Opt__Mag_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefMag_IntScheme",   HOFFSET(InputPara_t,Opt__RefMag_IntScheme  ), H5T_NATIVE_INT              );
#  endif
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__Pot_IntScheme",      HOFFSET(InputPara_t,Opt__Pot_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Rho_IntScheme",      HOFFSET(InputPara_t,Opt__Rho_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Gra_IntScheme",      HOFFSET(InputPara_t,Opt__Gra_IntScheme     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RefPot_IntScheme",   HOFFSET(InputPara_t,Opt__RefPot_IntScheme  ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "IntMonoCoeff",            HOFFSET(InputPara_t,IntMonoCoeff           ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "IntOppSign0thOrder",      HOFFSET(InputPara_t,IntOppSign0thOrder     ), H5T_NATIVE_INT              );

// data dump
   H5Tinsert( H5_TypeID, "Opt__Output_Total",       HOFFSET(InputPara_t,Opt__Output_Total      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Part",        HOFFSET(InputPara_t,Opt__Output_Part       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_User",        HOFFSET(InputPara_t,Opt__Output_User       ), H5T_NATIVE_INT              );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Output_ParText",     HOFFSET(InputPara_t,Opt__Output_ParText    ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Output_BasePS",      HOFFSET(InputPara_t,Opt__Output_BasePS     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Base",        HOFFSET(InputPara_t,Opt__Output_Base       ), H5T_NATIVE_INT              );
#  ifdef GRAVITY
   H5Tinsert( H5_TypeID, "Opt__Output_Pot",         HOFFSET(InputPara_t,Opt__Output_Pot        ), H5T_NATIVE_INT              );
#  endif
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Output_ParDens",     HOFFSET(InputPara_t,Opt__Output_ParDens    ), H5T_NATIVE_INT              );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Output_CC_Mag",      HOFFSET(InputPara_t,Opt__Output_CC_Mag     ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Output_Mode",        HOFFSET(InputPara_t,Opt__Output_Mode       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Step",        HOFFSET(InputPara_t,Opt__Output_Step       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Output_Dt",          HOFFSET(InputPara_t,Opt__Output_Dt         ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Output_PartX",            HOFFSET(InputPara_t,Output_PartX           ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Output_PartY",            HOFFSET(InputPara_t,Output_PartY           ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "Output_PartZ",            HOFFSET(InputPara_t,Output_PartZ           ), H5T_NATIVE_DOUBLE           );
   H5Tinsert( H5_TypeID, "InitDumpID",              HOFFSET(InputPara_t,InitDumpID             ), H5T_NATIVE_INT              );

// miscellaneous
   H5Tinsert( H5_TypeID, "Opt__Verbose",            HOFFSET(InputPara_t,Opt__Verbose           ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingBarrier",      HOFFSET(InputPara_t,Opt__TimingBarrier     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingBalance",      HOFFSET(InputPara_t,Opt__TimingBalance     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__TimingMPI",          HOFFSET(InputPara_t,Opt__TimingMPI         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordNote",         HOFFSET(InputPara_t,Opt__RecordNote        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordUnphy",        HOFFSET(InputPara_t,Opt__RecordUnphy       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordMemory",       HOFFSET(InputPara_t,Opt__RecordMemory      ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordPerformance",  HOFFSET(InputPara_t,Opt__RecordPerformance ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__ManualControl",      HOFFSET(InputPara_t,Opt__ManualControl     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__RecordUser",         HOFFSET(InputPara_t,Opt__RecordUser        ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__OptimizeAggressive", HOFFSET(InputPara_t,Opt__OptimizeAggressive), H5T_NATIVE_INT              );

// simulation checks
   H5Tinsert( H5_TypeID, "Opt__Ck_Refine",          HOFFSET(InputPara_t,Opt__Ck_Refine         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_ProperNesting",   HOFFSET(InputPara_t,Opt__Ck_ProperNesting  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Conservation",    HOFFSET(InputPara_t,Opt__Ck_Conservation   ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_NormPassive",     HOFFSET(InputPara_t,Opt__Ck_NormPassive    ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Restrict",        HOFFSET(InputPara_t,Opt__Ck_Restrict       ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_Finite",          HOFFSET(InputPara_t,Opt__Ck_Finite         ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_PatchAllocate",   HOFFSET(InputPara_t,Opt__Ck_PatchAllocate  ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_FluxAllocate",    HOFFSET(InputPara_t,Opt__Ck_FluxAllocate   ), H5T_NATIVE_INT              );
#  if ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "Opt__Ck_Negative",        HOFFSET(InputPara_t,Opt__Ck_Negative       ), H5T_NATIVE_INT              );
#  endif
   H5Tinsert( H5_TypeID, "Opt__Ck_MemFree",         HOFFSET(InputPara_t,Opt__Ck_MemFree        ), H5T_NATIVE_DOUBLE           );
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "Opt__Ck_Particle",        HOFFSET(InputPara_t,Opt__Ck_Particle       ), H5T_NATIVE_INT              );
#  endif
#  ifdef MHD
   H5Tinsert( H5_TypeID, "Opt__Ck_InterfaceB",      HOFFSET(InputPara_t,Opt__Ck_InterfaceB     ), H5T_NATIVE_INT              );
   H5Tinsert( H5_TypeID, "Opt__Ck_DivergenceB",     HOFFSET(InputPara_t,Opt__Ck_DivergenceB    ), H5T_NATIVE_INT              );
#  endif

// flag tables
#  if ( NLEVEL > 1 )
   H5Tinsert( H5_TypeID, "FlagTable_Rho",          HOFFSET(InputPara_t,FlagTable_Rho           ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_RhoGradient",  HOFFSET(InputPara_t,FlagTable_RhoGradient   ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_Lohner",       HOFFSET(InputPara_t,FlagTable_Lohner        ), H5_TypeID_Arr_NLvM1_4Double );

// store the user-defined thresholds at all levels
   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
//    key for each level
      sprintf( Key, "FlagTable_User_Lv%02d", lv );

//    assuming the offset between successive FlagTable_User pointers is "PtrSize_hvl", which is equal to "sizeof( hvl_t )"
      H5Tinsert( H5_TypeID, Key, HOFFSET(InputPara_t,FlagTable_User)+lv*PtrSize_hvl, H5_TypeID_Arr_NLvM1_VLDouble );
   }

#  if   ( MODEL == HYDRO )
   H5Tinsert( H5_TypeID, "FlagTable_PresGradient", HOFFSET(InputPara_t,FlagTable_PresGradient  ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_EngyGradient", HOFFSET(InputPara_t,FlagTable_EngyGradient  ), H5_TypeID_Arr_NLvM1Double   );         
   H5Tinsert( H5_TypeID, "FlagTable_Vorticity",    HOFFSET(InputPara_t,FlagTable_Vorticity     ), H5_TypeID_Arr_NLvM1Double   );
   H5Tinsert( H5_TypeID, "FlagTable_Jeans",        HOFFSET(InputPara_t,FlagTable_Jeans         ), H5_TypeID_Arr_NLvM1Double   );
#  ifdef MHD
   H5Tinsert( H5_TypeID, "FlagTable_Current",      HOFFSET(InputPara_t,FlagTable_Current       ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  ifdef SRHD
   H5Tinsert( H5_TypeID, "FlagTable_LorentzFactorGradient", HOFFSET(InputPara_t,FlagTable_LorentzFactorGradient  ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  elif ( MODEL == ELBDM )
   H5Tinsert( H5_TypeID, "FlagTable_EngyDensity",  HOFFSET(InputPara_t,FlagTable_EngyDensity   ), H5_TypeID_Arr_NLvM1_2Double );
#  endif
#  ifdef PARTICLE
   H5Tinsert( H5_TypeID, "FlagTable_NParPatch",    HOFFSET(InputPara_t,FlagTable_NParPatch     ), H5_TypeID_Arr_NLvM1Int      );
   H5Tinsert( H5_TypeID, "FlagTable_NParCell",     HOFFSET(InputPara_t,FlagTable_NParCell      ), H5_TypeID_Arr_NLvM1Int      );
   H5Tinsert( H5_TypeID, "FlagTable_ParMassCell",  HOFFSET(InputPara_t,FlagTable_ParMassCell   ), H5_TypeID_Arr_NLvM1Double   );
#  endif
#  endif


// free memory
   H5_Status = H5Tclose( H5_TypeID_Arr_3Int           );
   H5_Status = H5Tclose( H5_TypeID_Arr_6Int           );
#  if ( NCOMP_PASSIVE > 0 )
   H5_Status = H5Tclose( H5_TypeID_Arr_NPassive       );
#  endif
#  if ( NLEVEL > 1 )
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1Int       );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1Double    );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_2Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_4Double  );
   H5_Status = H5Tclose( H5_TypeID_Arr_NLvM1_VLDouble );
#  endif
   H5_Status = H5Tclose( H5_TypeID_Arr_3Double        );
   H5_Status = H5Tclose( H5_TypeID_VarStr             );

} // FUNCTION : GetCompound_InputPara



#endif // #ifdef SUPPORT_HDF5
