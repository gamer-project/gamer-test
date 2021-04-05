#ifndef __GLOBAL_H__
#define __GLOBAL_H__



#include "Macro.h"
#include "AMR.h"
#ifdef SUPPORT_LIBYT
#include <libyt.h>
#endif


// **********************************************************************************************************
// ** Variables in CAPITAL letters are loaded from the parameter file "Input__Parameter". Please refer to  **
// ** that file for detailed descritpions of each variable.                                                **
// **********************************************************************************************************


// 1. common global variables
// ============================================================================================================
extern AMR_t     *amr;                                // AMR tree structure
extern double     Time[NLEVEL];                       // "present"  physical time at each level
extern double     Time_Prev[NLEVEL];                  // "previous" physical time at each level
extern double     dTime_AllLv[NLEVEL];                // current evolution physical time interval at each level
extern long       AdvanceCounter[NLEVEL];             // number of sub-steps that each level has been evolved
extern long       NCorrUnphy[NLEVEL];                 // number of cells corrected by either OPT__1ST_FLUX_CORR or MIN_DENS/PRES
extern long       Step;                               // number of main steps
extern double     dTime_Base;                         // physical time interval at the base level

extern double     FlagTable_Rho        [NLEVEL-1];    // refinement criterion of density
extern double     FlagTable_RhoGradient[NLEVEL-1];    // refinement criterion of density gradient
extern double     FlagTable_Lohner     [NLEVEL-1][4]; // refinement criterion based on Lohner's error estimator
extern double    *FlagTable_User       [NLEVEL-1];    // user-defined refinement criterion
extern double    *DumpTable;                          // dump table recording the physical times to output data
extern int        DumpTable_NDump;                    // number of data dumps in the dump table
extern int        DumpID;                             // index of the current output file
extern double     DumpTime;                           // time of the next dump (for OPT__OUTPUT_MODE=1)

extern int        MPI_Rank;                           // MPI rank ID in the MPI_COMM_WORLD
extern int        MPI_Rank_X[3];                      // order of this MPI process in x/y/z directions
extern int        MPI_SibRank[26];                    // sibling MPI rank ID: same order as sibling[26]
extern int        NX0[3];                             //  of base-level cells per process in x/y/z directions
extern int        NPatchTotal[NLEVEL];                // total number of patches in all ranks
extern int       *BaseP;                              // table recording the IDs of the base-level patches
extern int        Flu_ParaBuf;                        // number of parallel buffers to exchange all fluid
                                                      // variables for the fluid solver and fluid refinement

extern int        PassiveNorm_NVar, PassiveNorm_VarIdx[NCOMP_PASSIVE];
extern int        PassiveIntFrac_NVar, PassiveIntFrac_VarIdx[NCOMP_PASSIVE];

extern double     BOX_SIZE, DT__MAX, DT__FLUID, DT__FLUID_INIT, END_T, OUTPUT_DT, DT__SYNC_PARENT_LV, DT__SYNC_CHILDREN_LV;
extern long int   END_STEP;
extern int        NX0_TOT[3], OUTPUT_STEP, REGRID_COUNT, FLU_GPU_NPGROUP, SRC_GPU_NPGROUP, OMP_NTHREAD;
extern int        MPI_NRank, MPI_NRank_X[3];
extern int        GPU_NSTREAM, FLAG_BUFFER_SIZE, FLAG_BUFFER_SIZE_MAXM1_LV, FLAG_BUFFER_SIZE_MAXM2_LV, MAX_LEVEL;

extern int        OPT__UM_IC_LEVEL, OPT__UM_IC_NVAR, OPT__UM_IC_LOAD_NRANK, OPT__GPUID_SELECT, OPT__PATCH_COUNT;
extern int        INIT_DUMPID, INIT_SUBSAMPLING_NCELL, OPT__TIMING_BARRIER, OPT__REUSE_MEMORY, RESTART_LOAD_NRANK;
extern double     OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z, AUTO_REDUCE_DT_FACTOR, AUTO_REDUCE_DT_FACTOR_MIN;
extern double     OPT__CK_MEMFREE, INT_MONO_COEFF, UNIT_L, UNIT_M, UNIT_T, UNIT_V, UNIT_D, UNIT_E, UNIT_P;
extern bool       OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_USER, OPT__FLAG_LOHNER_DENS, OPT__FLAG_REGION;
extern int        OPT__FLAG_USER_NUM;
extern bool       OPT__DT_USER, OPT__RECORD_DT, OPT__RECORD_MEMORY, OPT__MEMORY_POOL, OPT__RESTART_RESET;
extern bool       OPT__FIXUP_RESTRICT, OPT__INIT_RESTRICT, OPT__VERBOSE, OPT__MANUAL_CONTROL, OPT__UNIT;
extern bool       OPT__INT_TIME, OPT__OUTPUT_USER, OPT__OUTPUT_BASE, OPT__OVERLAP_MPI, OPT__TIMING_BALANCE;
extern bool       OPT__OUTPUT_BASEPS, OPT__CK_REFINE, OPT__CK_PROPER_NESTING, OPT__CK_FINITE, OPT__RECORD_PERFORMANCE;
extern bool       OPT__CK_RESTRICT, OPT__CK_PATCH_ALLOCATE, OPT__FIXUP_FLUX, OPT__CK_FLUX_ALLOCATE, OPT__CK_NORMALIZE_PASSIVE;
extern bool       OPT__UM_IC_DOWNGRADE, OPT__UM_IC_REFINE, OPT__TIMING_MPI;
extern bool       OPT__CK_CONSERVATION, OPT__RESET_FLUID, OPT__RECORD_USER, OPT__NORMALIZE_PASSIVE, AUTO_REDUCE_DT;
extern bool       OPT__OPTIMIZE_AGGRESSIVE, OPT__INIT_GRID_WITH_OMP, OPT__NO_FLAG_NEAR_BOUNDARY;
extern bool       OPT__RECORD_NOTE, OPT__RECORD_UNPHY, INT_OPP_SIGN_0TH_ORDER;
extern bool       OPT__INT_FRAC_PASSIVE_LR;

extern UM_IC_Format_t     OPT__UM_IC_FORMAT;
extern TestProbID_t       TESTPROB_ID;
extern OptInit_t          OPT__INIT;
extern IntScheme_t        OPT__FLU_INT_SCHEME, OPT__REF_FLU_INT_SCHEME;
extern OptOutputFormat_t  OPT__OUTPUT_TOTAL;
extern OptOutputPart_t    OPT__OUTPUT_PART;
extern OptOutputMode_t    OPT__OUTPUT_MODE;
extern OptFluBC_t         OPT__BC_FLU[6];          // boundary conditions of fluid at (-x,+x,-y,+y,-z,+z) faces
extern OptLohnerForm_t    OPT__FLAG_LOHNER_FORM;
extern OptCorrAfterSync_t OPT__CORR_AFTER_ALL_SYNC;
extern OptTimeStepLevel_t OPT__DT_LEVEL;



// 2. global variables for different applications
// ============================================================================================================
// (2-1) fluid solver in different models
#if   ( MODEL == HYDRO )
extern double           FlagTable_PresGradient[NLEVEL-1], FlagTable_Vorticity[NLEVEL-1], FlagTable_Jeans[NLEVEL-1];
extern double           GAMMA, MINMOD_COEFF, MOLECULAR_WEIGHT, ISO_TEMP;
extern LR_Limiter_t     OPT__LR_LIMITER;
extern Opt1stFluxCorr_t OPT__1ST_FLUX_CORR;
extern OptRSolver1st_t  OPT__1ST_FLUX_CORR_SCHEME;
extern bool             OPT__FLAG_PRES_GRADIENT, OPT__FLAG_LOHNER_ENGY, OPT__FLAG_LOHNER_PRES, OPT__FLAG_LOHNER_TEMP;
extern bool             OPT__FLAG_VORTICITY, OPT__FLAG_JEANS, JEANS_MIN_PRES, OPT__LAST_RESORT_FLOOR;
extern bool             OPT__OUTPUT_DIVVEL, OPT__OUTPUT_MACH, OPT__OUTPUT_PRES, OPT__OUTPUT_CS;
extern bool             OPT__OUTPUT_TEMP;
extern int              OPT__CK_NEGATIVE, JEANS_MIN_PRES_LEVEL, JEANS_MIN_PRES_NCELL;
extern double           MIN_DENS, MIN_PRES, MIN_EINT, MIN_TEMP;
#ifdef DUAL_ENERGY
extern double           DUAL_ENERGY_SWITCH;
#endif
#ifdef MHD
extern double           FlagTable_Current[NLEVEL-1];
extern IntScheme_t      OPT__MAG_INT_SCHEME, OPT__REF_MAG_INT_SCHEME;
extern bool             OPT__FIXUP_ELECTRIC, OPT__CK_INTERFACE_B, OPT__OUTPUT_CC_MAG, OPT__FLAG_CURRENT;
extern bool             OPT__OUTPUT_DIVMAG;
extern int              OPT__CK_DIVERGENCE_B;
extern double           UNIT_B;
extern bool             OPT__INIT_BFIELD_BYFILE;
#endif

#elif ( MODEL == ELBDM )
extern double           DT__PHASE, FlagTable_EngyDensity[NLEVEL-1][2];
extern bool             OPT__FLAG_ENGY_DENSITY, OPT__INT_PHASE, ELBDM_TAYLOR3_AUTO;
extern double           ELBDM_TAYLOR3_COEFF, ELBDM_MASS, ELBDM_PLANCK_CONST, ELBDM_ETA, MIN_DENS;
#ifdef QUARTIC_SELF_INTERACTION
extern double           ELBDM_LAMBDA;
#endif

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// (2-2) self-gravity
// ============================================================================================================
#ifdef GRAVITY
extern double        AveDensity_Init;     // initial average mass density (in all levels)
extern int           Pot_ParaBuf;         // number of parallel buffers to exchange potential for the
                                          // Poisson/Gravity solvers and the potential refinement
extern int           Rho_ParaBuf;         // number of parallel buffers to exchange density for the Poisson solver
extern real         *GreenFuncK;
extern double        GFUNC_COEFF0;
extern double        DT__GRAVITY;
extern double        NEWTON_G;
extern int           POT_GPU_NPGROUP;
extern bool          OPT__OUTPUT_POT, OPT__GRA_P5_GRADIENT, OPT__SELF_GRAVITY, OPT__GRAVITY_EXTRA_MASS;
extern double        SOR_OMEGA;
extern int           SOR_MAX_ITER, SOR_MIN_ITER;
extern double        MG_TOLERATED_ERROR;
extern int           MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH;
extern char          EXT_POT_TABLE_NAME[MAX_STRING];
extern double        EXT_POT_TABLE_DH, EXT_POT_TABLE_EDGEL[3];
extern int           EXT_POT_TABLE_NPOINT[3], EXT_POT_TABLE_FLOAT8;
extern IntScheme_t   OPT__POT_INT_SCHEME, OPT__RHO_INT_SCHEME, OPT__GRA_INT_SCHEME, OPT__REF_POT_INT_SCHEME;
extern OptPotBC_t    OPT__BC_POT;
extern OptExtAcc_t   OPT__EXT_ACC;
extern OptExtPot_t   OPT__EXT_POT;

extern double ExtAcc_AuxArray    [EXT_ACC_NAUX_MAX];
extern double ExtPot_AuxArray_Flt[EXT_POT_NAUX_MAX];
extern int    ExtPot_AuxArray_Int[EXT_POT_NAUX_MAX];
extern ExtAcc_t CPUExtAcc_Ptr;
extern ExtPot_t CPUExtPot_Ptr;
#ifdef GPU
extern ExtAcc_t GPUExtAcc_Ptr;
extern ExtPot_t GPUExtPot_Ptr;
#endif
#endif // #ifdef GRAVITY


// (2-3) cosmology simulations
// ============================================================================================================
#ifdef COMOVING
extern double     A_INIT, OMEGA_M0, DT__MAX_DELTA_A, HUBBLE0;
#endif


// (2-4) load balance
// ============================================================================================================
#ifdef LOAD_BALANCE
extern double     LB_INPUT__WLI_MAX;                  // LB->WLI_Max loaded from "Input__Parameter"
#ifdef PARTICLE
extern double     LB_INPUT__PAR_WEIGHT;               // LB->Par_Weight loaded from "Input__Parameter"
#endif
extern bool       OPT__RECORD_LOAD_BALANCE;
#endif
extern bool       OPT__MINIMIZE_MPI_BARRIER;


// (2-5) particle
// ============================================================================================================
#ifdef PARTICLE
extern double          DT__PARVEL, DT__PARVEL_MAX, DT__PARACC;
extern bool            OPT__OUTPUT_PAR_TEXT, OPT__CK_PARTICLE, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL;
extern int             OPT__PARTICLE_COUNT, OPT__FLAG_NPAR_PATCH, FlagTable_NParPatch[NLEVEL-1], FlagTable_NParCell[NLEVEL-1];
extern double          FlagTable_ParMassCell[NLEVEL-1];
extern ParOutputDens_t OPT__OUTPUT_PAR_DENS;
#endif


// (2-6) yt inline analysis
// ============================================================================================================
#ifdef SUPPORT_LIBYT
extern char            YT_SCRIPT[MAX_STRING];
extern yt_verbose      YT_VERBOSE;
#endif


// (2-7) Grackle
// ============================================================================================================
#ifdef SUPPORT_GRACKLE
extern bool            GRACKLE_ACTIVATE;
extern bool            GRACKLE_VERBOSE;
extern bool            GRACKLE_COOLING;
extern GracklePriChe_t GRACKLE_PRIMORDIAL;
extern bool            GRACKLE_METAL;
extern bool            GRACKLE_UV;
extern bool            GRACKLE_CMB_FLOOR;
extern bool            GRACKLE_PE_HEATING;
extern double          GRACKLE_PE_HEATING_RATE;
extern char            GRACKLE_CLOUDY_TABLE[MAX_STRING];
extern int             GRACKLE_THREE_BODY_RATE;
extern bool            GRACKLE_CIE_COOLING;
extern int             GRACKLE_H2_OPA_APPROX;
extern int             CHE_GPU_NPGROUP;
#endif


// (2-8) star formation
// ============================================================================================================
#ifdef STAR_FORMATION
extern SF_CreateStarScheme_t SF_CREATE_STAR_SCHEME;
extern int                   SF_CREATE_STAR_RSEED;
extern int                   SF_CREATE_STAR_DET_RANDOM;
extern int                   SF_CREATE_STAR_MIN_LEVEL;
extern double                SF_CREATE_STAR_MIN_GAS_DENS;
extern double                SF_CREATE_STAR_MASS_EFF;
extern double                SF_CREATE_STAR_MIN_STAR_MASS;
extern double                SF_CREATE_STAR_MAX_STAR_MFRAC;
#endif


// (2-9) equation of state
// =======================================================================================================
#if ( MODEL == HYDRO )
extern double EoS_AuxArray_Flt[EOS_NAUX_MAX];
extern int    EoS_AuxArray_Int[EOS_NAUX_MAX];
extern EoS_DE2P_t EoS_DensEint2Pres_CPUPtr;
extern EoS_DP2E_t EoS_DensPres2Eint_CPUPtr;
extern EoS_DP2C_t EoS_DensPres2CSqr_CPUPtr;
extern EoS_DE2T_t EoS_DensEint2Temp_CPUPtr;
extern EoS_DT2P_t EoS_DensTemp2Pres_CPUPtr;
extern EoS_GENE_t EoS_General_CPUPtr;
#ifdef GPU
extern EoS_DE2P_t EoS_DensEint2Pres_GPUPtr;
extern EoS_DP2E_t EoS_DensPres2Eint_GPUPtr;
extern EoS_DP2C_t EoS_DensPres2CSqr_GPUPtr;
extern EoS_DE2T_t EoS_DensEint2Temp_GPUPtr;
extern EoS_DT2P_t EoS_DensTemp2Pres_GPUPtr;
extern EoS_GENE_t EoS_General_GPUPtr;
#endif
extern EoS_t EoS;
#endif // HYDRO


// (2-10) source terms
// =======================================================================================================
extern SrcTerms_t SrcTerms;
#if ( MODEL == HYDRO )
extern double     Src_Dlep_AuxArray_Flt[SRC_NAUX_DLEP];
extern int        Src_Dlep_AuxArray_Int[SRC_NAUX_DLEP];
#endif
extern double     Src_User_AuxArray_Flt[SRC_NAUX_USER];
extern int        Src_User_AuxArray_Int[SRC_NAUX_USER];


// (2-11) user-defined derived fields
// =======================================================================================================
extern bool OPT__OUTPUT_USER_FIELD;
extern int  UserDerField_Num;
extern char (*UserDerField_Label)[MAX_STRING];
extern char (*UserDerField_Unit )[MAX_STRING];
// Flu_DerivedField_User_Ptr is defined in Flu_DerivedField_User.cpp
extern void (*Flu_DerivedField_User_Ptr)( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                                          const int NCellInX, const int NCellInY, const int NCellInZ,
                                          const int NGhost, const double dh );


// (2-12) feedback
// =======================================================================================================
#ifdef FEEDBACK
extern int  FB_LEVEL, FB_RSEED;
extern bool FB_SNE, FB_USER;
extern bool FB_Any;
#endif



// 3. CPU (host) arrays for transferring data between CPU and GPU
// ============================================================================================================
extern real       (*h_Flu_Array_F_In [2])[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real       (*h_Flu_Array_F_Out[2])[FLU_NOUT][ CUBE(PS2) ];
extern real       (*h_Flux_Array[2])[9][NFLUX_TOTAL][ SQR(PS2) ];
extern double     (*h_Corner_Array_F [2])[3];
#ifdef DUAL_ENERGY
extern char       (*h_DE_Array_F_Out [2])[ CUBE(PS2) ];
#endif
#ifdef MHD
extern real       (*h_Mag_Array_F_In [2])[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real       (*h_Mag_Array_F_Out[2])[NCOMP_MAG][ PS2P1*SQR(PS2)          ];
extern real       (*h_Ele_Array      [2])[9][NCOMP_ELE][ PS2P1*PS2 ];
#endif

#ifdef GRAVITY
extern real       (*h_Rho_Array_P     [2])[RHO_NXT][RHO_NXT][RHO_NXT];
extern real       (*h_Pot_Array_P_In  [2])[POT_NXT][POT_NXT][POT_NXT];
extern real       (*h_Pot_Array_P_Out [2])[GRA_NXT][GRA_NXT][GRA_NXT];
extern real       (*h_Flu_Array_G     [2])[GRA_NIN][PS1][PS1][PS1];
extern double     (*h_Corner_Array_PGT[2])[3];
#ifdef DUAL_ENERGY
extern char       (*h_DE_Array_G      [2])[PS1][PS1][PS1];
#endif
#ifdef MHD
extern real       (*h_Emag_Array_G    [2])[PS1][PS1][PS1];
#endif
extern real        *h_ExtPotTable;
extern void       **h_ExtPotGenePtr;

#ifdef UNSPLIT_GRAVITY
extern real       (*h_Pot_Array_USG_F [2])[ CUBE(USG_NXT_F) ];
extern real       (*h_Pot_Array_USG_G [2])[USG_NXT_G ][USG_NXT_G ][USG_NXT_G ];
extern real       (*h_Flu_Array_USG_G [2])[GRA_NIN-1][PS1][PS1][PS1];
#endif
#endif // #ifdef GRAVITY

#ifdef SUPPORT_GRACKLE
extern real       (*h_Che_Array[2]);
// do not declare Grackle variables for CUDA source files since they do not include <grackle.h>
#ifndef __CUDACC__
extern grackle_field_data *Che_FieldData;
extern code_units Che_Units;
#endif
#endif

extern real        *h_dt_Array_T[2];
extern real       (*h_Flu_Array_T[2])[FLU_NIN_T][ CUBE(PS1) ];
#ifdef GRAVITY
extern real       (*h_Pot_Array_T[2])[ CUBE(GRA_NXT) ];
#endif
#ifdef MHD
extern real       (*h_Mag_Array_T[2])[NCOMP_MAG][ PS1P1*SQR(PS1) ];
#endif

#if ( MODEL == HYDRO )
extern real        *h_EoS_Table[EOS_NTABLE_MAX];
#endif

extern real       (*h_Flu_Array_S_In [2])[FLU_NIN_S ][ CUBE(SRC_NXT)           ];
extern real       (*h_Flu_Array_S_Out[2])[FLU_NOUT_S][ CUBE(PS1)               ];
#ifdef MHD
extern real       (*h_Mag_Array_S_In [2])[NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ];
#endif
extern double     (*h_Corner_Array_S[2])[3];
#if ( MODEL == HYDRO )
extern real       (*h_SrcDlepProf_Data)[SRC_DLEP_PROF_NBINMAX];
extern real        *h_SrcDlepProf_Radius;
#endif



// 4/5. GPU (device) global memory arrays and timers
// ============================================================================================================
/*** These global variables are NOT included here. Instead, they are included by individual files
     only if necessary. ***/



// 6. global variables related to different fields
// ============================================================================================================
/*** Defined in Field.h ***/



#endif // #ifndef __GLOBAL_H__
