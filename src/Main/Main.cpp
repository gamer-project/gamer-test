// define DEFINE_GLOBAL since this file **defines** all global variables
#define DEFINE_GLOBAL
#include "GAMER.h"
#undef DEFINE_GLOBAL
#include "CUFLU.h"

#ifdef GRAVITY
#include "CUPOT.h"
#endif



// ***********************
// **  GLOBAL VARIABLES **
// ***********************

// 1. common global variables
// =======================================================================================================
AMR_t               *amr = NULL;

double               Time[NLEVEL]           = { 0.0 };
double               dTime_AllLv[NLEVEL]    = { 0.0 };
long                 AdvanceCounter[NLEVEL] = { 0 };
long                 NCorrUnphy[NLEVEL]     = { 0 };
long                 Step                   = 0;
int                  DumpID                 = 0;
double               DumpTime               = 0.0;

double               dTime_Base;
double               Time_Prev            [NLEVEL];
double               FlagTable_Rho        [NLEVEL-1];
double               FlagTable_RhoGradient[NLEVEL-1];
double               FlagTable_Lohner     [NLEVEL-1][4];
double              *FlagTable_User       [NLEVEL-1];
double              *DumpTable = NULL;
int                  DumpTable_NDump;
int                  PassiveNorm_NVar, PassiveNorm_VarIdx[NCOMP_PASSIVE];
int                  PassiveIntFrac_NVar, PassiveIntFrac_VarIdx[NCOMP_PASSIVE];

int                  MPI_Rank, MPI_Rank_X[3], MPI_SibRank[26], NX0[3], NPatchTotal[NLEVEL];
int                 *BaseP = NULL;
int                  Flu_ParaBuf;

double               BOX_SIZE, DT__MAX, DT__FLUID, DT__FLUID_INIT, END_T, OUTPUT_DT, DT__SYNC_PARENT_LV, DT__SYNC_CHILDREN_LV;
long                 END_STEP;
int                  NX0_TOT[3], OUTPUT_STEP, REGRID_COUNT, FLU_GPU_NPGROUP, SRC_GPU_NPGROUP, OMP_NTHREAD;
int                  MPI_NRank, MPI_NRank_X[3];
int                  GPU_NSTREAM, FLAG_BUFFER_SIZE, FLAG_BUFFER_SIZE_MAXM1_LV, FLAG_BUFFER_SIZE_MAXM2_LV, MAX_LEVEL;

IntScheme_t          OPT__FLU_INT_SCHEME, OPT__REF_FLU_INT_SCHEME;
double               OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z, AUTO_REDUCE_DT_FACTOR, AUTO_REDUCE_DT_FACTOR_MIN;
double               OPT__CK_MEMFREE, INT_MONO_COEFF, UNIT_L, UNIT_M, UNIT_T, UNIT_V, UNIT_D, UNIT_E, UNIT_P;
int                  OPT__UM_IC_LEVEL, OPT__UM_IC_NVAR, OPT__UM_IC_LOAD_NRANK, OPT__GPUID_SELECT, OPT__PATCH_COUNT;
int                  INIT_DUMPID, INIT_SUBSAMPLING_NCELL, OPT__TIMING_BARRIER, OPT__REUSE_MEMORY, RESTART_LOAD_NRANK;
bool                 OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_USER, OPT__FLAG_LOHNER_DENS, OPT__FLAG_REGION;
int                  OPT__FLAG_USER_NUM;
bool                 OPT__DT_USER, OPT__RECORD_DT, OPT__RECORD_MEMORY, OPT__MEMORY_POOL, OPT__RESTART_RESET;
bool                 OPT__FIXUP_RESTRICT, OPT__INIT_RESTRICT, OPT__VERBOSE, OPT__MANUAL_CONTROL, OPT__UNIT;
bool                 OPT__INT_TIME, OPT__OUTPUT_USER, OPT__OUTPUT_BASE, OPT__OVERLAP_MPI, OPT__TIMING_BALANCE;
bool                 OPT__OUTPUT_BASEPS, OPT__CK_REFINE, OPT__CK_PROPER_NESTING, OPT__CK_FINITE, OPT__RECORD_PERFORMANCE;
bool                 OPT__CK_RESTRICT, OPT__CK_PATCH_ALLOCATE, OPT__FIXUP_FLUX, OPT__CK_FLUX_ALLOCATE, OPT__CK_NORMALIZE_PASSIVE;
bool                 OPT__UM_IC_DOWNGRADE, OPT__UM_IC_REFINE, OPT__TIMING_MPI;
bool                 OPT__CK_CONSERVATION, OPT__RESET_FLUID, OPT__RECORD_USER, OPT__NORMALIZE_PASSIVE, AUTO_REDUCE_DT;
bool                 OPT__OPTIMIZE_AGGRESSIVE, OPT__INIT_GRID_WITH_OMP, OPT__NO_FLAG_NEAR_BOUNDARY;
bool                 OPT__RECORD_NOTE, OPT__RECORD_UNPHY, INT_OPP_SIGN_0TH_ORDER;
bool                 OPT__INT_FRAC_PASSIVE_LR;

UM_IC_Format_t       OPT__UM_IC_FORMAT;
TestProbID_t         TESTPROB_ID;
OptInit_t            OPT__INIT;
OptOutputFormat_t    OPT__OUTPUT_TOTAL;
OptOutputPart_t      OPT__OUTPUT_PART;
OptOutputMode_t      OPT__OUTPUT_MODE;
OptFluBC_t           OPT__BC_FLU[6];
OptLohnerForm_t      OPT__FLAG_LOHNER_FORM;
OptCorrAfterSync_t   OPT__CORR_AFTER_ALL_SYNC;
OptTimeStepLevel_t   OPT__DT_LEVEL;


// 2. global variables for different applications
// =======================================================================================================
// (2-1) fluid solver in different models
#if   ( MODEL == HYDRO )
double               FlagTable_PresGradient[NLEVEL-1], FlagTable_Vorticity[NLEVEL-1], FlagTable_Jeans[NLEVEL-1];
double               GAMMA, MINMOD_COEFF, MOLECULAR_WEIGHT, ISO_TEMP;
LR_Limiter_t         OPT__LR_LIMITER;
Opt1stFluxCorr_t     OPT__1ST_FLUX_CORR;
OptRSolver1st_t      OPT__1ST_FLUX_CORR_SCHEME;
bool                 OPT__FLAG_PRES_GRADIENT, OPT__FLAG_LOHNER_ENGY, OPT__FLAG_LOHNER_PRES, OPT__FLAG_LOHNER_TEMP;
bool                 OPT__FLAG_VORTICITY, OPT__FLAG_JEANS, JEANS_MIN_PRES, OPT__LAST_RESORT_FLOOR;
bool                 OPT__OUTPUT_DIVVEL, OPT__OUTPUT_MACH, OPT__OUTPUT_PRES, OPT__OUTPUT_CS;
bool                 OPT__OUTPUT_TEMP;
int                  OPT__CK_NEGATIVE, JEANS_MIN_PRES_LEVEL, JEANS_MIN_PRES_NCELL;
double               MIN_DENS, MIN_PRES, MIN_EINT, MIN_TEMP;
#ifdef DUAL_ENERGY
double               DUAL_ENERGY_SWITCH;
#endif
#ifdef MHD
double               FlagTable_Current[NLEVEL-1];
IntScheme_t          OPT__MAG_INT_SCHEME, OPT__REF_MAG_INT_SCHEME;
bool                 OPT__FIXUP_ELECTRIC, OPT__CK_INTERFACE_B, OPT__OUTPUT_CC_MAG, OPT__FLAG_CURRENT;
bool                 OPT__OUTPUT_DIVMAG;
int                  OPT__CK_DIVERGENCE_B;
double               UNIT_B;
bool                 OPT__INIT_BFIELD_BYFILE;
#endif

#elif ( MODEL == ELBDM )
double               DT__PHASE, FlagTable_EngyDensity[NLEVEL-1][2];
bool                 OPT__FLAG_ENGY_DENSITY, OPT__INT_PHASE;
bool                 ELBDM_TAYLOR3_AUTO;
double               ELBDM_TAYLOR3_COEFF;
double               ELBDM_MASS, ELBDM_PLANCK_CONST, ELBDM_ETA, MIN_DENS;
#ifdef QUARTIC_SELF_INTERACTION
double               ELBDM_LAMBDA;
#endif

#else
#error : unsupported MODEL !!
#endif // MODEL

// (2-2) self-gravity
#ifdef GRAVITY
double               AveDensity_Init = -1.0;    // initialize it as <= 0 to check if it is properly set later
int                  Pot_ParaBuf, Rho_ParaBuf;

real                *GreenFuncK = NULL;
double               GFUNC_COEFF0;
double               DT__GRAVITY;
double               NEWTON_G;
int                  POT_GPU_NPGROUP;
bool                 OPT__OUTPUT_POT, OPT__GRA_P5_GRADIENT, OPT__SELF_GRAVITY, OPT__GRAVITY_EXTRA_MASS;
double               SOR_OMEGA;
int                  SOR_MAX_ITER, SOR_MIN_ITER;
double               MG_TOLERATED_ERROR;
int                  MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH;
char                 EXT_POT_TABLE_NAME[MAX_STRING];
double               EXT_POT_TABLE_DH, EXT_POT_TABLE_EDGEL[3];
int                  EXT_POT_TABLE_NPOINT[3], EXT_POT_TABLE_FLOAT8;
IntScheme_t          OPT__POT_INT_SCHEME, OPT__RHO_INT_SCHEME, OPT__GRA_INT_SCHEME, OPT__REF_POT_INT_SCHEME;
OptPotBC_t           OPT__BC_POT;
OptExtAcc_t          OPT__EXT_ACC;
OptExtPot_t          OPT__EXT_POT;

// external gravity variables
// a. auxiliary arrays
double ExtAcc_AuxArray    [EXT_ACC_NAUX_MAX];
double ExtPot_AuxArray_Flt[EXT_POT_NAUX_MAX];
int    ExtPot_AuxArray_Int[EXT_POT_NAUX_MAX];

// b. function pointers
ExtAcc_t CPUExtAcc_Ptr = NULL;
ExtPot_t CPUExtPot_Ptr = NULL;
#ifdef GPU
ExtAcc_t GPUExtAcc_Ptr = NULL;
ExtPot_t GPUExtPot_Ptr = NULL;
#endif
#endif // #ifdef GRAVITY

// (2-3) cosmological simulations
#ifdef COMOVING
double               A_INIT, OMEGA_M0, DT__MAX_DELTA_A, HUBBLE0;
#endif

// (2-4) load balance
#ifdef LOAD_BALANCE
double               LB_INPUT__WLI_MAX;
#ifdef PARTICLE
double               LB_INPUT__PAR_WEIGHT;
#endif
bool                 OPT__RECORD_LOAD_BALANCE;
#endif
bool                 OPT__MINIMIZE_MPI_BARRIER;

// (2-5) particle
#ifdef PARTICLE
double               DT__PARVEL, DT__PARVEL_MAX, DT__PARACC;
bool                 OPT__OUTPUT_PAR_TEXT, OPT__CK_PARTICLE, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL;
int                  OPT__PARTICLE_COUNT, OPT__FLAG_NPAR_PATCH, FlagTable_NParPatch[NLEVEL-1], FlagTable_NParCell[NLEVEL-1];
double               FlagTable_ParMassCell[NLEVEL-1];
ParOutputDens_t      OPT__OUTPUT_PAR_DENS;
#endif

// (2-6) yt inline analysis
#ifdef SUPPORT_LIBYT
char                 YT_SCRIPT[MAX_STRING];
yt_verbose           YT_VERBOSE;
#endif

// (2-7) Grackle
#ifdef SUPPORT_GRACKLE
bool                 GRACKLE_ACTIVATE;
bool                 GRACKLE_VERBOSE;
bool                 GRACKLE_COOLING;
GracklePriChe_t      GRACKLE_PRIMORDIAL;
bool                 GRACKLE_METAL;
bool                 GRACKLE_UV;
bool                 GRACKLE_CMB_FLOOR;
bool                 GRACKLE_PE_HEATING;
double               GRACKLE_PE_HEATING_RATE;
char                 GRACKLE_CLOUDY_TABLE[MAX_STRING];
int                  GRACKLE_THREE_BODY_RATE;
bool                 GRACKLE_CIE_COOLING;
int                  GRACKLE_H2_OPA_APPROX;
int                  CHE_GPU_NPGROUP;
#endif

// (2-8) star formation
#ifdef STAR_FORMATION
SF_CreateStarScheme_t SF_CREATE_STAR_SCHEME;
int                   SF_CREATE_STAR_RSEED;
int                   SF_CREATE_STAR_DET_RANDOM;
int                   SF_CREATE_STAR_MIN_LEVEL;
double                SF_CREATE_STAR_MIN_GAS_DENS;
double                SF_CREATE_STAR_MASS_EFF;
double                SF_CREATE_STAR_MIN_STAR_MASS;
double                SF_CREATE_STAR_MAX_STAR_MFRAC;
#endif

// (2-9) equation of state
#if ( MODEL == HYDRO )
// a. auxiliary arrays
double EoS_AuxArray_Flt[EOS_NAUX_MAX];
int    EoS_AuxArray_Int[EOS_NAUX_MAX];

// b. function pointers
EoS_DE2P_t EoS_DensEint2Pres_CPUPtr = NULL;
EoS_DP2E_t EoS_DensPres2Eint_CPUPtr = NULL;
EoS_DP2C_t EoS_DensPres2CSqr_CPUPtr = NULL;
EoS_DE2T_t EoS_DensEint2Temp_CPUPtr = NULL;
EoS_DT2P_t EoS_DensTemp2Pres_CPUPtr = NULL;
EoS_GENE_t EoS_General_CPUPtr       = NULL;
#ifdef GPU
EoS_DE2P_t EoS_DensEint2Pres_GPUPtr = NULL;
EoS_DP2E_t EoS_DensPres2Eint_GPUPtr = NULL;
EoS_DP2C_t EoS_DensPres2CSqr_GPUPtr = NULL;
EoS_DE2T_t EoS_DensEint2Temp_GPUPtr = NULL;
EoS_DT2P_t EoS_DensTemp2Pres_GPUPtr = NULL;
EoS_GENE_t EoS_General_GPUPtr       = NULL;
#endif

// c. data structure for the CPU/GPU solvers
EoS_t EoS;
#endif // HYDRO

// (2-10) source terms
SrcTerms_t SrcTerms;
#if ( MODEL == HYDRO )
double     Src_Dlep_AuxArray_Flt[SRC_NAUX_DLEP];
int        Src_Dlep_AuxArray_Int[SRC_NAUX_DLEP];
#endif
double     Src_User_AuxArray_Flt[SRC_NAUX_USER];
int        Src_User_AuxArray_Int[SRC_NAUX_USER];

// (2-11) user-defined derived fields
bool OPT__OUTPUT_USER_FIELD;
int  UserDerField_Num                  = -1;    // must be negative for Output_DumpData_Total_HDF5()
char (*UserDerField_Label)[MAX_STRING] = NULL;
char (*UserDerField_Unit )[MAX_STRING] = NULL;

// (2-12) feedback
#ifdef FEEDBACK
int  FB_LEVEL, FB_RSEED;
bool FB_SNE, FB_USER;
bool FB_Any;
#endif


// 3. CPU (host) arrays for transferring data between CPU and GPU
// =======================================================================================================
// (3-1) fluid solver
real (*h_Flu_Array_F_In [2])[FLU_NIN ][ CUBE(FLU_NXT) ]            = { NULL, NULL };
real (*h_Flu_Array_F_Out[2])[FLU_NOUT][ CUBE(PS2) ]                = { NULL, NULL };
real (*h_Flux_Array[2])[9][NFLUX_TOTAL][ SQR(PS2) ]                = { NULL, NULL };
double (*h_Corner_Array_F[2])[3]                                   = { NULL, NULL };
#ifdef DUAL_ENERGY
char (*h_DE_Array_F_Out[2])[ CUBE(PS2) ]                           = { NULL, NULL };
#endif
#ifdef MHD
real (*h_Mag_Array_F_In [2])[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = { NULL, NULL };
real (*h_Mag_Array_F_Out[2])[NCOMP_MAG][ PS2P1*SQR(PS2)          ] = { NULL, NULL };
real (*h_Ele_Array      [2])[9][NCOMP_ELE][ PS2P1*PS2 ]            = { NULL, NULL };
#endif
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
real (*h_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ]  = NULL;
real (*h_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ]  = NULL;
real (*h_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ]  = NULL;
real (*h_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ]  = NULL;
#ifdef MHD
real (*h_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ]        = NULL;
real (*h_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ]        = NULL;
#endif
#endif // FLU_SCHEME

#ifdef GRAVITY
// (3-2) Poisson and gravity solver
real   (*h_Rho_Array_P     [2])[RHO_NXT][RHO_NXT][RHO_NXT]         = { NULL, NULL };
real   (*h_Pot_Array_P_In  [2])[POT_NXT][POT_NXT][POT_NXT]         = { NULL, NULL };
real   (*h_Pot_Array_P_Out [2])[GRA_NXT][GRA_NXT][GRA_NXT]         = { NULL, NULL };
real   (*h_Flu_Array_G     [2])[GRA_NIN][PS1][PS1][PS1]            = { NULL, NULL };
double (*h_Corner_Array_PGT[2])[3]                                 = { NULL, NULL };
#ifdef DUAL_ENERGY
char   (*h_DE_Array_G      [2])[PS1][PS1][PS1]                     = { NULL, NULL };
#endif
#ifdef MHD
real   (*h_Emag_Array_G    [2])[PS1][PS1][PS1]                     = { NULL, NULL };
#endif
real    *h_ExtPotTable                                             = NULL;
void   **h_ExtPotGenePtr                                           = NULL;

// (3-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*h_Pot_Array_USG_F[2])[ CUBE(USG_NXT_F) ]                    = { NULL, NULL };
real (*h_Pot_Array_USG_G[2])[USG_NXT_G][USG_NXT_G][USG_NXT_G]      = { NULL, NULL };
real (*h_Flu_Array_USG_G[2])[GRA_NIN-1][PS1][PS1][PS1]             = { NULL, NULL };
#endif
#endif // #ifdef GRAVITY

// (3-4) Grackle chemistry
#ifdef SUPPORT_GRACKLE
real (*h_Che_Array[2])                                             = { NULL, NULL };
grackle_field_data *Che_FieldData                                  = NULL;
code_units Che_Units;
#endif

// (3-5) dt solver
real  *h_dt_Array_T[2]                                             = { NULL, NULL };
real (*h_Flu_Array_T[2])[FLU_NIN_T][ CUBE(PS1) ]                   = { NULL, NULL };
#ifdef GRAVITY
real (*h_Pot_Array_T[2])[ CUBE(GRA_NXT) ]                          = { NULL, NULL };
#endif
#ifdef MHD
real (*h_Mag_Array_T[2])[NCOMP_MAG][ PS1P1*SQR(PS1) ]              = { NULL, NULL };
#endif

// (3-6) EoS tables
#if ( MODEL == HYDRO )
real *h_EoS_Table[EOS_NTABLE_MAX];
#endif

// (3-7) source terms
real (*h_Flu_Array_S_In [2])[FLU_NIN_S ][ CUBE(SRC_NXT)  ]         = { NULL, NULL };
real (*h_Flu_Array_S_Out[2])[FLU_NOUT_S][ CUBE(PS1)      ]         = { NULL, NULL };
#ifdef MHD
real (*h_Mag_Array_S_In [2])[NCOMP_MAG][ SRC_NXT_P1*SQR(SRC_NXT) ] = { NULL, NULL };
#endif
double (*h_Corner_Array_S[2])[3]                                   = { NULL, NULL };
#if ( MODEL == HYDRO )
real (*h_SrcDlepProf_Data)[SRC_DLEP_PROF_NBINMAX]                  = NULL;
real  *h_SrcDlepProf_Radius                                        = NULL;
#endif


// 4. GPU (device) global memory arrays
// =======================================================================================================
#ifdef GPU
// (4-1) fluid solver
real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ]               = NULL;
real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ]                   = NULL;
real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ]                   = NULL;
double (*d_Corner_Array_F)[3]                                      = NULL;
#ifdef DUAL_ENERGY
char (*d_DE_Array_F_Out)[ PS2*PS2*PS2 ]                            = NULL;
#endif
#ifdef MHD
real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ]    = NULL;
real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)          ]    = NULL;
real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ]               = NULL;
#endif
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
real (*d_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ]  = NULL;
real (*d_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ]  = NULL;
real (*d_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ]  = NULL;
real (*d_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ]  = NULL;
#ifdef MHD
real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ]        = NULL;
real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ]        = NULL;
#endif
#endif // FLU_SCHEME

#ifdef GRAVITY
// (4-2) Poisson and gravity solver
real   (*d_Rho_Array_P    )[ CUBE(RHO_NXT) ]                       = NULL;
real   (*d_Pot_Array_P_In )[ CUBE(POT_NXT) ]                       = NULL;
real   (*d_Pot_Array_P_Out)[ CUBE(GRA_NXT) ]                       = NULL;
real   (*d_Flu_Array_G    )[GRA_NIN][ CUBE(PS1) ]                  = NULL;
double (*d_Corner_Array_PGT)[3]                                    = NULL;
#ifdef DUAL_ENERGY
char   (*d_DE_Array_G     )[ CUBE(PS1) ]                           = NULL;
#endif
#ifdef MHD
real   (*d_Emag_Array_G   )[ CUBE(PS1) ]                           = NULL;
#endif
real    *d_ExtPotTable                                             = NULL;
void   **d_ExtPotGenePtr                                           = NULL;

// (4-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ]                       = NULL;
real (*d_Pot_Array_USG_G)[ CUBE(USG_NXT_G) ]                       = NULL;
real (*d_Flu_Array_USG_G)[GRA_NIN-1][ CUBE(PS1) ]                  = NULL;
#endif
#endif // #ifdef GRAVITY

// (4-4) Grackle chemistry

// (4-5) dt solver
real *d_dt_Array_T                                                 = NULL;
real (*d_Flu_Array_T)[FLU_NIN_T][ CUBE(PS1) ]                      = NULL;
#ifdef GRAVITY
real (*d_Pot_Array_T)[ CUBE(GRA_NXT) ]                             = NULL;
#endif
#ifdef MHD
real (*d_Mag_Array_T)[NCOMP_MAG][ PS1P1*SQR(PS1) ]                 = NULL;
#endif

// (4-6) EoS tables
#if ( MODEL == HYDRO )
real *d_EoS_Table[EOS_NTABLE_MAX];
#endif

// (4-7) source terms
real (*d_Flu_Array_S_In )[FLU_NIN_S ][ CUBE(SRC_NXT)  ]            = NULL;
real (*d_Flu_Array_S_Out)[FLU_NOUT_S][ CUBE(PS1)      ]            = NULL;
#ifdef MHD
real (*d_Mag_Array_S_In)[NCOMP_MAG  ][ SRC_NXT_P1*SQR(SRC_NXT) ]   = NULL;
#endif
double (*d_Corner_Array_S)[3]                                      = NULL;
#if ( MODEL == HYDRO )
real (*d_SrcDlepProf_Data)[SRC_DLEP_PROF_NBINMAX]                  = NULL;
real  *d_SrcDlepProf_Radius                                        = NULL;
#endif

#endif // #ifdef GPU


// 5. timers
// =======================================================================================================
#ifdef TIMING
Timer_t *Timer_Main[7];
Timer_t *Timer_MPI[3];
Timer_t *Timer_dt         [NLEVEL];
Timer_t *Timer_Flu_Advance[NLEVEL];
Timer_t *Timer_Gra_Advance[NLEVEL];
Timer_t *Timer_Src_Advance[NLEVEL];
Timer_t *Timer_Che_Advance[NLEVEL];
Timer_t *Timer_SF         [NLEVEL];
Timer_t *Timer_FB_Advance [NLEVEL];
Timer_t *Timer_FixUp      [NLEVEL];
Timer_t *Timer_Flag       [NLEVEL];
Timer_t *Timer_Refine     [NLEVEL];
Timer_t *Timer_GetBuf     [NLEVEL][9];
Timer_t *Timer_Lv         [NLEVEL];
Timer_t *Timer_Par_Update [NLEVEL][3];
Timer_t *Timer_Par_2Sib   [NLEVEL];
Timer_t *Timer_Par_2Son   [NLEVEL];
Timer_t *Timer_Par_Collect[NLEVEL];
Timer_t *Timer_Par_MPI    [NLEVEL][6];
#endif

#ifdef TIMING_SOLVER
Timer_t *Timer_Pre         [NLEVEL][NSOLVER];
Timer_t *Timer_Sol         [NLEVEL][NSOLVER];
Timer_t *Timer_Clo         [NLEVEL][NSOLVER];
Timer_t *Timer_Poi_PreRho  [NLEVEL];
Timer_t *Timer_Poi_PreFlu  [NLEVEL];
Timer_t *Timer_Poi_PrePot_C[NLEVEL];
Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif


// function pointer for recording the user-specified info
extern void (*Aux_Record_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  GAMER main function
//-------------------------------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{

// initialization
// ======================================================================================================
   Timer_t Timer_Total;
   Timer_Total.Start();

#  ifdef TIMING
   Timer_t  Timer_Init, Timer_Other;
   Timer_Init.Start();
#  endif


   Init_GAMER( &argc, &argv );

   if ( OPT__RECORD_NOTE )
   {
      Aux_TakeNote();

#     ifdef GPU
      CUAPI_DiagnoseDevice();
#     endif
   }

   Output_DumpData( 0 );

   if ( OPT__PATCH_COUNT > 0 )            Aux_Record_PatchCount();
   if ( OPT__RECORD_MEMORY )              Aux_GetMemInfo();
   if ( OPT__RECORD_USER ) {
      if ( Aux_Record_User_Ptr != NULL )  Aux_Record_User_Ptr();
      else
         Aux_Error( ERROR_INFO, "Aux_Record_User_Ptr == NULL for OPT__RECORD_USER !!\n" );
   }

#  ifdef PARTICLE
   if ( OPT__PARTICLE_COUNT > 0 )         Par_Aux_Record_ParticleCount();
#  endif

   Aux_Check();

#  ifdef TIMING
   Aux_ResetTimer();
#  endif


#  ifdef TIMING
   Timer_Init.Stop();
#  endif
// ======================================================================================================



// main loop
// ======================================================================================================
   MPI_Barrier( MPI_COMM_WORLD );


   while ( (Time[0]-END_T < -1.e-10)  &&  (Step < END_STEP) )
   {

#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Main[0]->Start();    // timer for one iteration
#     endif


//    1. advance all physical attributes by one global time-step
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   EvolveLevel( 0, NULL_REAL ),     Timer_Main[2],   TIMER_ON   );

      Step ++;
//    ---------------------------------------------------------------------------------------------------


//    2. apply various corrections
//       --> synchronize particles, restrict data, recalculate potential and particle acceleration, ...
//    ---------------------------------------------------------------------------------------------------
      if ( OPT__CORR_AFTER_ALL_SYNC == CORR_AFTER_SYNC_EVERY_STEP )
      TIMING_FUNC(   Flu_CorrAfterAllSync(),          Timer_Main[6],   TIMER_ON   );
//    ---------------------------------------------------------------------------------------------------


//    3. output data and execute auxiliary functions
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   Output_DumpData( 1 ),            Timer_Main[3],   TIMER_ON   );

      if ( OPT__PATCH_COUNT == 1 )
      TIMING_FUNC(   Aux_Record_PatchCount(),         Timer_Main[4],   TIMER_ON   );

      if ( OPT__RECORD_MEMORY )
      TIMING_FUNC(   Aux_GetMemInfo(),                Timer_Main[4],   TIMER_ON   );

      if ( OPT__RECORD_USER )
      TIMING_FUNC(   Aux_Record_User_Ptr(),           Timer_Main[4],   TIMER_ON   );

      if ( OPT__RECORD_UNPHY )
      TIMING_FUNC(   Aux_Record_CorrUnphy(),          Timer_Main[4],   TIMER_ON   );

#     ifdef PARTICLE
      if ( OPT__PARTICLE_COUNT == 1 )
      TIMING_FUNC(   Par_Aux_Record_ParticleCount(),  Timer_Main[4],   TIMER_ON   );
#     endif

      TIMING_FUNC(   Aux_Check(),                     Timer_Main[4],   TIMER_ON   );
//    ---------------------------------------------------------------------------------------------------


//    4. perform yt inline analysis
//    ---------------------------------------------------------------------------------------------------
#     ifdef SUPPORT_LIBYT
      YT_Inline();
#     endif
//    ---------------------------------------------------------------------------------------------------


//    5. check whether to manually terminate the run
//    ---------------------------------------------------------------------------------------------------
      int Terminate = false;

//    enable this functionality only if OPT__MANUAL_CONTROL is on
      if ( OPT__MANUAL_CONTROL )
      TIMING_FUNC(   End_StopManually( Terminate ),   Timer_Main[4],   TIMER_ON   );
//    ---------------------------------------------------------------------------------------------------


//    6. check whether to redistribute all patches for LOAD_BALANCE
//    ---------------------------------------------------------------------------------------------------
#     ifdef LOAD_BALANCE
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );
#     ifdef TIMING
      Timer_Main[5]->Start();    // timer for load balance
#     endif

      if ( LB_EstimateLoadImbalance() > amr->LB->WLI_Max )
      {
         if ( MPI_Rank == 0 )
         {
            Aux_Message( stdout, "Weighted load-imbalance factor (%13.7e) > threshold (%13.7e) ",
                         amr->LB->WLI, amr->LB->WLI_Max );
            Aux_Message( stdout, "--> redistributing all patches ...\n" );
         }

         const bool   Redistribute_Yes = true;
         const bool   ResetLB_Yes      = true;
#        ifdef PARTICLE
         const double ParWeight        = amr->LB->Par_Weight;
#        else
         const double ParWeight        = 0.0;
#        endif
         const int    AllLv            = -1;

         LB_Init_LoadBalance( Redistribute_Yes, ParWeight, ResetLB_Yes, AllLv );

         if ( OPT__PATCH_COUNT > 0 )         Aux_Record_PatchCount();

#        ifdef PARTICLE
         if ( OPT__PARTICLE_COUNT > 0 )      Par_Aux_Record_ParticleCount();
#        endif
      } // if ( LB_EstimateLoadImbalance() > amr->LB->WLI_Max )

#     ifdef TIMING
      Timer_Main[5]->Stop();
#     endif
#     endif // #ifdef LOAD_BALANCE
//    ---------------------------------------------------------------------------------------------------


//    7. record timing
//    ---------------------------------------------------------------------------------------------------
#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Main[0]->Stop();

      Timer_Other.Start();

      if ( OPT__RECORD_PERFORMANCE )
      Aux_Record_Performance( Timer_Main[0]->GetValue() );

      Aux_Record_Timing();

      Aux_ResetTimer();

      Timer_Other.Stop();
#     endif
//    ---------------------------------------------------------------------------------------------------


      if ( Terminate )  break;

   } // while ( (Time[0]-END_T < -1.e-10)  &&  (Step < END_STEP) )


   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Total.Stop();
// ======================================================================================================



// termination
// ======================================================================================================
// output the final result
   Output_DumpData( 2 );


// record the total simulation time
#  ifdef TIMING
   Aux_AccumulatedTiming( Timer_Total.GetValue(), Timer_Init.GetValue(), Timer_Other.GetValue() );
#  endif

   if ( MPI_Rank == 0  &&  OPT__RECORD_NOTE )
   {
      FILE *Note = fopen( "Record__Note", "a" );
      fprintf( Note, "\n" );
      fprintf( Note, "Total Processing Time : %lf s\n", Timer_Total.GetValue() );
      fprintf( Note, "\n" );
      fclose( Note );
   }


   End_GAMER();
// ======================================================================================================

   return 0;

} // FUNCTION : Main

