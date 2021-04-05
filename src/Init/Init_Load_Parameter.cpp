#include "GAMER.h"
#include "ReadPara.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_Parameter
// Description :  Load the runtime parameters
//
// Note        :  Currently the filename is fixed to "Input__Parameter"
//-------------------------------------------------------------------------------------------------------
void Init_Load_Parameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char  FileName[] = "Input__Parameter";
   ReadPara_t *ReadPara   = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_INPUT_TABLE",     &VARIABLE,                       DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************

// simulation scale
   ReadPara->Add( "BOX_SIZE",                   &BOX_SIZE,                       -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "NX0_TOT_X",                  &NX0_TOT[0],                     -1,               PS2,           NoMax_int      );
   ReadPara->Add( "NX0_TOT_Y",                  &NX0_TOT[1],                     -1,               PS2,           NoMax_int      );
   ReadPara->Add( "NX0_TOT_Z",                  &NX0_TOT[2],                     -1,               PS2,           NoMax_int      );
// do not check MPI_NRANK_X/Y/Z since they are deprecated
   ReadPara->Add( "MPI_NRANK_X",                &MPI_NRank_X[0],                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MPI_NRANK_Y",                &MPI_NRank_X[1],                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MPI_NRANK_Z",                &MPI_NRank_X[2],                 -1,               NoMin_int,     NoMax_int      );
// do not check OMP_NTHREAD since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OMP_NTHREAD",                &OMP_NTHREAD,                    -1,               NoMin_int,     NoMax_int      );
// do not check END_T and END_STEP since they may be reset by test problems or restart
   ReadPara->Add( "END_T",                      &END_T,                          -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "END_STEP",                   &END_STEP,                       -1L,              NoMin_long,    NoMax_long     );


// test problems
   ReadPara->Add( "TESTPROB_ID",                &TESTPROB_ID,                     0,               0,             NoMax_int      );


// code units
   ReadPara->Add( "OPT__UNIT",                  &OPT__UNIT,                       false,           Useless_bool,  Useless_bool   );
// do not check units since they are validated by Init_Unit()
   ReadPara->Add( "UNIT_L",                     &UNIT_L,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_M",                     &UNIT_M,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_T",                     &UNIT_T,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_V",                     &UNIT_V,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_D",                     &UNIT_D,                         -1.0,             NoMin_double,  NoMax_double   );


// boundary conditions
   ReadPara->Add( "OPT__BC_FLU_XM",             &OPT__BC_FLU[0],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_XP",             &OPT__BC_FLU[1],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_YM",             &OPT__BC_FLU[2],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_YP",             &OPT__BC_FLU[3],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_ZM",             &OPT__BC_FLU[4],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_ZP",             &OPT__BC_FLU[5],                 -1,               1,             4              );
#  ifdef GRAVITY
   ReadPara->Add( "OPT__BC_POT",                &OPT__BC_POT,                    -1,               1,             2              );
// do not check GFUNC_COEFF0 since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "GFUNC_COEFF0",               &GFUNC_COEFF0,                   -1.0,             NoMin_double,  NoMax_double   );
#  endif


// particle
#  ifdef PARTICLE
// do no check PAR_NPAR since it may be reset by restart
   ReadPara->Add( "PAR_NPAR",                   &amr->Par->NPar_Active_AllRank,  -1L,              NoMin_long,    NoMax_long     );
   ReadPara->Add( "PAR_INIT",                   &amr->Par->Init,                 -1,               1,             3              );
   ReadPara->Add( "PAR_IC_FORMAT",              &amr->Par->ParICFormat,      PAR_IC_FORMAT_ATT_ID, 1,             2              );
   ReadPara->Add( "PAR_IC_MASS",                &amr->Par->ParICMass,            -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "PAR_INTERP",                 &amr->Par->Interp,                PAR_INTERP_CIC,  1,             3              );
   ReadPara->Add( "PAR_INTEG",                  &amr->Par->Integ,                 PAR_INTEG_KDK,   1,             2              );
   ReadPara->Add( "PAR_IMPROVE_ACC",            &amr->Par->ImproveAcc,            true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "PAR_PREDICT_POS",            &amr->Par->PredictPos,            true,            Useless_bool,  Useless_bool   );
// do not check PAR_REMOVE_CELL since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "PAR_REMOVE_CELL",            &amr->Par->RemoveCell,           -1.0,             NoMin_double,  NoMax_double   );
#  endif // #ifdef PARTICLE


// cosmology
#  ifdef COMOVING
   ReadPara->Add( "A_INIT",                     &A_INIT,                         -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "OMEGA_M0",                   &OMEGA_M0,                       -1.0,             0.0,           1.0            );
   ReadPara->Add( "HUBBLE0",                    &HUBBLE0,                        -1.0,             Eps_double,    1.0            );
#  endif


// time-step
   ReadPara->Add( "DT__MAX",                    &DT__MAX,                        -1.0,             NoMin_double,  NoMax_double   );
// do not check DT__FLUID/FLUID_INIT/GRAVITY/PARVEL_MAX since they may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "DT__FLUID",                  &DT__FLUID,                      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "DT__FLUID_INIT",             &DT__FLUID_INIT,                 -1.0,             NoMin_double,  NoMax_double   );
#  ifdef GRAVITY
   ReadPara->Add( "DT__GRAVITY",                &DT__GRAVITY,                    -1.0,             NoMin_double,  NoMax_double   );
#  endif
#  if ( MODEL == ELBDM )
   ReadPara->Add( "DT__PHASE",                  &DT__PHASE,                       0.0,             0.0,           NoMax_double   );
#  endif
#  ifdef PARTICLE
   ReadPara->Add( "DT__PARVEL",                 &DT__PARVEL,                      0.5,             0.0,           NoMax_double   );
   ReadPara->Add( "DT__PARVEL_MAX",             &DT__PARVEL_MAX,                 -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "DT__PARACC",                 &DT__PARACC,                      0.5,             0.0,           NoMax_double   );
#  endif
#  ifdef COMOVING
   ReadPara->Add( "DT__MAX_DELTA_A",            &DT__MAX_DELTA_A,                 0.01,            0.0,           NoMax_double   );
#  endif
   ReadPara->Add( "DT__SYNC_PARENT_LV",         &DT__SYNC_PARENT_LV,              0.1,             0.0,           NoMax_double   );
   ReadPara->Add( "DT__SYNC_CHILDREN_LV",       &DT__SYNC_CHILDREN_LV,            0.1,             0.0,           1.0            );
   ReadPara->Add( "OPT__DT_USER",               &OPT__DT_USER,                    false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__DT_LEVEL",              &OPT__DT_LEVEL,                   3,               1,             3              );
   ReadPara->Add( "OPT__RECORD_DT",             &OPT__RECORD_DT,                  true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "AUTO_REDUCE_DT",             &AUTO_REDUCE_DT,                  true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "AUTO_REDUCE_DT_FACTOR",      &AUTO_REDUCE_DT_FACTOR,           0.8,             Eps_double,    1.0            );
   ReadPara->Add( "AUTO_REDUCE_DT_FACTOR_MIN",  &AUTO_REDUCE_DT_FACTOR_MIN,       0.1,             0.0,           1.0            );


// grid refinement
   ReadPara->Add( "REGRID_COUNT",               &REGRID_COUNT,                    4,               1,             NoMax_int      );
   ReadPara->Add( "FLAG_BUFFER_SIZE",           &FLAG_BUFFER_SIZE,                PS1,             0,             PS1            );
   ReadPara->Add( "FLAG_BUFFER_SIZE_MAXM1_LV",  &FLAG_BUFFER_SIZE_MAXM1_LV,      -1,               NoMin_int,     PS1            );
   ReadPara->Add( "FLAG_BUFFER_SIZE_MAXM2_LV",  &FLAG_BUFFER_SIZE_MAXM2_LV,      -1,               NoMin_int,     PS1            );
   ReadPara->Add( "MAX_LEVEL",                  &MAX_LEVEL,                       TOP_LEVEL,       0,             TOP_LEVEL      );
   ReadPara->Add( "OPT__FLAG_RHO",              &OPT__FLAG_RHO,                   false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_RHO_GRADIENT",     &OPT__FLAG_RHO_GRADIENT,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_DENS",      &OPT__FLAG_LOHNER_DENS,           false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO )
   ReadPara->Add( "OPT__FLAG_PRES_GRADIENT",    &OPT__FLAG_PRES_GRADIENT,         false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_VORTICITY",        &OPT__FLAG_VORTICITY,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_JEANS",            &OPT__FLAG_JEANS,                 false,           Useless_bool,  Useless_bool   );
#  ifdef MHD
   ReadPara->Add( "OPT__FLAG_CURRENT",          &OPT__FLAG_CURRENT,               false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__FLAG_LOHNER_ENGY",      &OPT__FLAG_LOHNER_ENGY,           false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_PRES",      &OPT__FLAG_LOHNER_PRES,           false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_TEMP",      &OPT__FLAG_LOHNER_TEMP,           false,           Useless_bool,  Useless_bool   );
#  endif
#  if ( MODEL == ELBDM )
   ReadPara->Add( "OPT__FLAG_ENGY_DENSITY",     &OPT__FLAG_ENGY_DENSITY,          false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__FLAG_LOHNER_FORM",      &OPT__FLAG_LOHNER_FORM,           LOHNER_FLASH2,   1,             4              );
   ReadPara->Add( "OPT__FLAG_USER",             &OPT__FLAG_USER,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_USER_NUM",         &OPT__FLAG_USER_NUM,              1,               1,             NoMax_int      );
   ReadPara->Add( "OPT__FLAG_REGION",           &OPT__FLAG_REGION,                false,           Useless_bool,  Useless_bool   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__FLAG_NPAR_PATCH",       &OPT__FLAG_NPAR_PATCH,            0,               0,             2              );
   ReadPara->Add( "OPT__FLAG_NPAR_CELL",        &OPT__FLAG_NPAR_CELL,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_PAR_MASS_CELL",    &OPT__FLAG_PAR_MASS_CELL,         false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__NO_FLAG_NEAR_BOUNDARY", &OPT__NO_FLAG_NEAR_BOUNDARY,      false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__PATCH_COUNT",           &OPT__PATCH_COUNT,                1,               0,             2              );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__PARTICLE_COUNT",        &OPT__PARTICLE_COUNT,             1,               0,             2              );
#  endif
   ReadPara->Add( "OPT__REUSE_MEMORY",          &OPT__REUSE_MEMORY,               2,               0,             2              );
   ReadPara->Add( "OPT__MEMORY_POOL",           &OPT__MEMORY_POOL,                false,           Useless_bool,  Useless_bool   );


// load balance
#  ifdef LOAD_BALANCE
   ReadPara->Add( "LB_INPUT__WLI_MAX",          &LB_INPUT__WLI_MAX,               0.1,             0.0,           NoMax_double   );
#  ifdef PARTICLE
   ReadPara->Add( "LB_INPUT__PAR_WEIGHT",       &LB_INPUT__PAR_WEIGHT,            0.0,             0.0,           NoMax_double   );
#  endif
   ReadPara->Add( "OPT__RECORD_LOAD_BALANCE",   &OPT__RECORD_LOAD_BALANCE,        true,            Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__MINIMIZE_MPI_BARRIER",  &OPT__MINIMIZE_MPI_BARRIER,       true,            Useless_bool,  Useless_bool   );


// source terms
   ReadPara->Add( "SRC_DELEPTONIZATION",        &SrcTerms.Deleptonization,        false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "SRC_USER",                   &SrcTerms.User,                   false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "SRC_GPU_NPGROUP",            &SRC_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );


// Grackle
#  ifdef SUPPORT_GRACKLE
   ReadPara->Add( "GRACKLE_ACTIVATE",           &GRACKLE_ACTIVATE,                true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_VERBOSE",            &GRACKLE_VERBOSE,                 true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_COOLING",            &GRACKLE_COOLING,                 true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_PRIMORDIAL",         &GRACKLE_PRIMORDIAL,              0,               0,             3              );
   ReadPara->Add( "GRACKLE_METAL",              &GRACKLE_METAL,                   false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_UV",                 &GRACKLE_UV,                      false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_CMB_FLOOR",          &GRACKLE_CMB_FLOOR,               true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_PE_HEATING",         &GRACKLE_PE_HEATING,              false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_PE_HEATING_RATE",    &GRACKLE_PE_HEATING_RATE,         8.5e-26,         0.0,           NoMax_double   );
   ReadPara->Add( "GRACKLE_CLOUDY_TABLE",        GRACKLE_CLOUDY_TABLE,            Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "GRACKLE_THREE_BODY_RATE",    &GRACKLE_THREE_BODY_RATE,         0,               0,             5              );
   ReadPara->Add( "GRACKLE_CIE_COOLING",        &GRACKLE_CIE_COOLING,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_H2_OPA_APPROX",      &GRACKLE_H2_OPA_APPROX,           0,               0,             1              );
// do not check CHE_GPU_NPGROUP since it may be reset by either Init_ResetDefaultParameter() or CUAPI_Set_Default_GPU_Parameter()
   ReadPara->Add( "CHE_GPU_NPGROUP",            &CHE_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );
#  endif


// star formation
#  ifdef STAR_FORMATION
   ReadPara->Add( "SF_CREATE_STAR_SCHEME",         &SF_CREATE_STAR_SCHEME,         0,              0,             1              );
   ReadPara->Add( "SF_CREATE_STAR_RSEED",          &SF_CREATE_STAR_RSEED,          123,            0,             NoMax_int      );
// do not check SF_CREATE_STAR_DET_RANDOM since its default depends on the makefile option BITWISE_REPRODUCIBILITY
   ReadPara->Add( "SF_CREATE_STAR_DET_RANDOM",     &SF_CREATE_STAR_DET_RANDOM,    -1,              NoMin_int,     NoMax_int      );
   ReadPara->Add( "SF_CREATE_STAR_MIN_LEVEL",      &SF_CREATE_STAR_MIN_LEVEL,      0,              NoMin_int,     TOP_LEVEL      );
   ReadPara->Add( "SF_CREATE_STAR_MIN_GAS_DENS",   &SF_CREATE_STAR_MIN_GAS_DENS,   1.0e1,          0.0,           NoMax_double   );
   ReadPara->Add( "SF_CREATE_STAR_MASS_EFF",       &SF_CREATE_STAR_MASS_EFF,       1.0e-2,         Eps_double,    1.0            );
   ReadPara->Add( "SF_CREATE_STAR_MIN_STAR_MASS",  &SF_CREATE_STAR_MIN_STAR_MASS,  0.0,            0.0,           NoMax_double   );
   ReadPara->Add( "SF_CREATE_STAR_MAX_STAR_MFRAC", &SF_CREATE_STAR_MAX_STAR_MFRAC, 0.5,            Eps_double,    1.0            );
#  endif


// feedback
#  ifdef FEEDBACK
   ReadPara->Add( "FB_LEVEL",                   &FB_LEVEL,                       -1,               NoMin_int,     TOP_LEVEL      );
   ReadPara->Add( "FB_RSEED",                   &FB_RSEED,                        456,             0,             NoMax_int      );
   ReadPara->Add( "FB_SNE",                     &FB_SNE,                          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "FB_USER",                    &FB_USER,                         false,           Useless_bool,  Useless_bool   );
#  endif


// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
#  if ( EOS == EOS_GAMMA )
   ReadPara->Add( "GAMMA",                      &GAMMA,                           5.0/3.0,         1.0,           NoMax_double   );
#  else
   ReadPara->Add( "GAMMA",                      &GAMMA,                           __DBL_MAX__,     NoMin_double,  NoMax_double   );
#  endif
   ReadPara->Add( "MOLECULAR_WEIGHT",           &MOLECULAR_WEIGHT,                0.6,             Eps_double,    NoMax_double   );
#  if ( EOS == EOS_ISOTHERMAL )
   ReadPara->Add( "ISO_TEMP",                   &ISO_TEMP,                       -1.0,             Eps_double,    NoMax_double   );
#  else
   ReadPara->Add( "ISO_TEMP",                   &ISO_TEMP,                       __DBL_MAX__,      NoMin_double,  NoMax_double   );
#  endif
   ReadPara->Add( "MINMOD_COEFF",               &MINMOD_COEFF,                    1.5,             1.0,           2.0            );
   ReadPara->Add( "OPT__LR_LIMITER",            &OPT__LR_LIMITER,                 VL_GMINMOD,      0,             5              );
   ReadPara->Add( "OPT__1ST_FLUX_CORR",         &OPT__1ST_FLUX_CORR,             -1,               NoMin_int,     2              );
#  ifdef MHD
   ReadPara->Add( "OPT__1ST_FLUX_CORR_SCHEME",  &OPT__1ST_FLUX_CORR_SCHEME,   RSOLVER_1ST_DEFAULT, NoMin_int,     4              );
#  else
   ReadPara->Add( "OPT__1ST_FLUX_CORR_SCHEME",  &OPT__1ST_FLUX_CORR_SCHEME,   RSOLVER_1ST_DEFAULT, NoMin_int,     3              );
#  endif
#  ifdef DUAL_ENERGY
   ReadPara->Add( "DUAL_ENERGY_SWITCH",         &DUAL_ENERGY_SWITCH,              2.0e-2,          0.0,           NoMax_double   );
#  endif
#  endif // #if ( MODEL == HYDRO )


// fluid solver in ELBDM
#  if ( MODEL == ELBDM )
   ReadPara->Add( "ELBDM_MASS",                 &ELBDM_MASS,                     -1.0,             Eps_double,    NoMax_double   );
// do not check ELBDM_PLANCK_CONST since it may be reset by Init_Unit()
   ReadPara->Add( "ELBDM_PLANCK_CONST",         &ELBDM_PLANCK_CONST,             -1.0,             NoMin_double,  NoMax_double   );
#  ifdef QUARTIC_SELF_INTERACTION
   ReadPara->Add( "ELBDM_LAMBDA",               &ELBDM_LAMBDA,                    1.0,             NoMin_double,  NoMax_double   );
#  endif
   ReadPara->Add( "ELBDM_TAYLOR3_COEFF",        &ELBDM_TAYLOR3_COEFF,             1.0/6.0,         NoMin_double,  NoMax_double   );
   ReadPara->Add( "ELBDM_TAYLOR3_AUTO",         &ELBDM_TAYLOR3_AUTO,              true,            Useless_bool,  Useless_bool   );
#  endif // #if ( MODEL == ELBDM )


// fluid solvers in all models
// do not check FLU_GPU_NPGROUP and GPU_NSTREAM since they may be reset by either Init_ResetDefaultParameter() or CUAPI_Set_Default_GPU_Parameter()
   ReadPara->Add( "FLU_GPU_NPGROUP",            &FLU_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "GPU_NSTREAM",                &GPU_NSTREAM,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__FIXUP_FLUX",            &OPT__FIXUP_FLUX,                 true,            Useless_bool,  Useless_bool   );
#  ifdef MHD
   ReadPara->Add( "OPT__FIXUP_ELECTRIC",        &OPT__FIXUP_ELECTRIC,             true,            Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__FIXUP_RESTRICT",        &OPT__FIXUP_RESTRICT,             true,            Useless_bool,  Useless_bool   );
// do not check OPT__CORR_AFTER_ALL_SYNC since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OPT__CORR_AFTER_ALL_SYNC",   &OPT__CORR_AFTER_ALL_SYNC,       -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__NORMALIZE_PASSIVE",     &OPT__NORMALIZE_PASSIVE,          true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__INT_FRAC_PASSIVE_LR",   &OPT__INT_FRAC_PASSIVE_LR,        true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OVERLAP_MPI",           &OPT__OVERLAP_MPI,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RESET_FLUID",           &OPT__RESET_FLUID,                false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO )
   ReadPara->Add( "OPT__LAST_RESORT_FLOOR",     &OPT__LAST_RESORT_FLOOR,          true,            Useless_bool,  Useless_bool   );
#  endif
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   ReadPara->Add( "MIN_DENS",                   &MIN_DENS,                        0.0,             0.0,           NoMax_double   );
#  endif
#  if ( MODEL == HYDRO )
   ReadPara->Add( "MIN_PRES",                   &MIN_PRES,                        0.0,             0.0,           NoMax_double   );
   ReadPara->Add( "MIN_EINT",                   &MIN_EINT,                        0.0,             0.0,           NoMax_double   );
   ReadPara->Add( "MIN_TEMP",                   &MIN_TEMP,                        0.0,             0.0,           NoMax_double   );
   ReadPara->Add( "JEANS_MIN_PRES",             &JEANS_MIN_PRES,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "JEANS_MIN_PRES_LEVEL",       &JEANS_MIN_PRES_LEVEL,           -1,               NoMin_int,     NLEVEL-1       );
   ReadPara->Add( "JEANS_MIN_PRES_NCELL",       &JEANS_MIN_PRES_NCELL,            4,               1,             NoMax_int      );
#  endif


// self-gravity
#  ifdef GRAVITY
// do not check NEWTON_G since it may be reset by Init_Unit()
   ReadPara->Add( "NEWTON_G",                   &NEWTON_G,                       -1.0,             NoMin_double,  NoMax_double   );
// do not check SOR_XXX since they may be reset by Init_Set_Default_SOR_Parameter()
   ReadPara->Add( "SOR_OMEGA",                  &SOR_OMEGA,                      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "SOR_MAX_ITER",               &SOR_MAX_ITER,                   -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "SOR_MIN_ITER",               &SOR_MIN_ITER,                   -1,               NoMin_int,     NoMax_int      );
// do not check MG_XXX since they may be reset by Init_Set_Default_MG_Parameter()
   ReadPara->Add( "MG_MAX_ITER",                &MG_MAX_ITER,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_NPRE_SMOOTH",             &MG_NPRE_SMOOTH,                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_NPOST_SMOOTH",            &MG_NPOST_SMOOTH,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_TOLERATED_ERROR",         &MG_TOLERATED_ERROR,             -1.0,             NoMin_double,  NoMax_double   );
// do not check POT_GPU_NPGROUP since it may be reset by either Init_ResetDefaultParameter() or CUAPI_Set_Default_GPU_Parameter()
   ReadPara->Add( "POT_GPU_NPGROUP",            &POT_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__GRA_P5_GRADIENT",       &OPT__GRA_P5_GRADIENT,            false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__SELF_GRAVITY",          &OPT__SELF_GRAVITY,               true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__EXT_ACC",               &OPT__EXT_ACC,                    0,               0,             1              );
   ReadPara->Add( "OPT__EXT_POT",               &OPT__EXT_POT,                    0,               0,             2              );
// do not check the parameters of external potential table here --> do it in Init_LoadExtPotTable()
   ReadPara->Add( "EXT_POT_TABLE_NAME",          EXT_POT_TABLE_NAME,              Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "EXT_POT_TABLE_NPOINT_X",     &EXT_POT_TABLE_NPOINT[0],        -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "EXT_POT_TABLE_NPOINT_Y",     &EXT_POT_TABLE_NPOINT[1],        -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "EXT_POT_TABLE_NPOINT_Z",     &EXT_POT_TABLE_NPOINT[2],        -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "EXT_POT_TABLE_DH",           &EXT_POT_TABLE_DH,               -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "EXT_POT_TABLE_EDGEL_X",      &EXT_POT_TABLE_EDGEL[0],          NoDef_double,    NoMin_double,  NoMax_double   );
   ReadPara->Add( "EXT_POT_TABLE_EDGEL_Y",      &EXT_POT_TABLE_EDGEL[1],          NoDef_double,    NoMin_double,  NoMax_double   );
   ReadPara->Add( "EXT_POT_TABLE_EDGEL_Z",      &EXT_POT_TABLE_EDGEL[2],          NoDef_double,    NoMin_double,  NoMax_double   );
// fix EXT_POT_TABLE_FLOAT8 to -1 for now since this option is not supported yet
   ReadPara->Add( "EXT_POT_TABLE_FLOAT8",       &EXT_POT_TABLE_FLOAT8,           -1,              -1,            -1              );
   ReadPara->Add( "OPT__GRAVITY_EXTRA_MASS",    &OPT__GRAVITY_EXTRA_MASS,         false,           Useless_bool,  Useless_bool   );
#  endif // #ifdef GRAVITY


// initialization
   ReadPara->Add( "OPT__INIT",                  &OPT__INIT,                      -1,               1,             3              );
   ReadPara->Add( "RESTART_LOAD_NRANK",         &RESTART_LOAD_NRANK,              1,               1,             NoMax_int      );
   ReadPara->Add( "OPT__RESTART_RESET",         &OPT__RESTART_RESET,              false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__UM_IC_LEVEL",           &OPT__UM_IC_LEVEL,                0,               0,             TOP_LEVEL      );
// do not check OPT__UM_IC_NVAR since it depends on OPT__INIT and MODEL
// --> also, we do not load the density field for ELBDM
#  if ( MODEL == ELBDM )
   ReadPara->Add( "OPT__UM_IC_NVAR",            &OPT__UM_IC_NVAR,                -1,               NoMin_int,     NCOMP_TOTAL-1  );
#  else
   ReadPara->Add( "OPT__UM_IC_NVAR",            &OPT__UM_IC_NVAR,                -1,               NoMin_int,     NCOMP_TOTAL    );
#  endif
   ReadPara->Add( "OPT__UM_IC_FORMAT",          &OPT__UM_IC_FORMAT,             UM_IC_FORMAT_VZYX, 1,             2              );
   ReadPara->Add( "OPT__UM_IC_DOWNGRADE",       &OPT__UM_IC_DOWNGRADE,            true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__UM_IC_REFINE",          &OPT__UM_IC_REFINE,               true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__UM_IC_LOAD_NRANK",      &OPT__UM_IC_LOAD_NRANK,           1,               1,             NoMax_int      );
   ReadPara->Add( "OPT__INIT_RESTRICT",         &OPT__INIT_RESTRICT,              true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__INIT_GRID_WITH_OMP",    &OPT__INIT_GRID_WITH_OMP,         true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__GPUID_SELECT",          &OPT__GPUID_SELECT,              -1,              -3,             NoMax_int      );
   ReadPara->Add( "INIT_SUBSAMPLING_NCELL",     &INIT_SUBSAMPLING_NCELL,          0,               0,             NoMax_int      );
#  ifdef MHD
   ReadPara->Add( "OPT__INIT_BFIELD_BYFILE",    &OPT__INIT_BFIELD_BYFILE,         false,           Useless_bool,  Useless_bool   );
#  endif

// interpolation schemes
   ReadPara->Add( "OPT__INT_TIME",              &OPT__INT_TIME,                   true,            Useless_bool,  Useless_bool   );
#  if ( MODEL == ELBDM )
   ReadPara->Add( "OPT__INT_PHASE",             &OPT__INT_PHASE,                  true,            Useless_bool,  Useless_bool   );
#  endif
// do not check OPT__FLU_INT_SCHEME and OPT__REF_FLU_INT_SCHEME since they may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OPT__FLU_INT_SCHEME",        &OPT__FLU_INT_SCHEME,             INT_DEFAULT,     NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__REF_FLU_INT_SCHEME",    &OPT__REF_FLU_INT_SCHEME,         INT_DEFAULT,     NoMin_int,     NoMax_int      );
#  ifdef MHD
   ReadPara->Add( "OPT__MAG_INT_SCHEME",        &OPT__MAG_INT_SCHEME,             INT_CQUAD,       NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__REF_MAG_INT_SCHEME",    &OPT__REF_MAG_INT_SCHEME,         INT_CQUAD,       NoMin_int,     NoMax_int      );
#  endif
#  ifdef GRAVITY
   ReadPara->Add( "OPT__POT_INT_SCHEME",        &OPT__POT_INT_SCHEME,             INT_QUAD,        4,             5              );
   ReadPara->Add( "OPT__RHO_INT_SCHEME",        &OPT__RHO_INT_SCHEME,             INT_CQUAD,       1,             7              );
   ReadPara->Add( "OPT__GRA_INT_SCHEME",        &OPT__GRA_INT_SCHEME,             INT_QUAD,        1,             7              );
   ReadPara->Add( "OPT__REF_POT_INT_SCHEME",    &OPT__REF_POT_INT_SCHEME,         INT_QUAD,        1,             7              );
#  endif
   ReadPara->Add( "INT_MONO_COEFF",             &INT_MONO_COEFF,                  2.0,             1.0,           4.0            );
#  if   ( MODEL == HYDRO )
   ReadPara->Add( "INT_OPP_SIGN_0TH_ORDER",     &INT_OPP_SIGN_0TH_ORDER,          true,            Useless_bool,  Useless_bool   );
#  elif ( MODEL == ELBDM )
   ReadPara->Add( "INT_OPP_SIGN_0TH_ORDER",     &INT_OPP_SIGN_0TH_ORDER,          false,           Useless_bool,  Useless_bool   );
#  else
#  error : unsupported MODEL !!
#  endif


// data dump
   ReadPara->Add( "OPT__OUTPUT_TOTAL",          &OPT__OUTPUT_TOTAL,               1,               0,             2              );
   ReadPara->Add( "OPT__OUTPUT_PART",           &OPT__OUTPUT_PART,                0,               0,             7              );
   ReadPara->Add( "OPT__OUTPUT_USER",           &OPT__OUTPUT_USER,                false,           Useless_bool,  Useless_bool   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__OUTPUT_PAR_TEXT",       &OPT__OUTPUT_PAR_TEXT,            false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__OUTPUT_BASEPS",         &OPT__OUTPUT_BASEPS,              false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_BASE",           &OPT__OUTPUT_BASE,                false,           Useless_bool,  Useless_bool   );
#  ifdef MHD
   ReadPara->Add( "OPT__OUTPUT_CC_MAG",         &OPT__OUTPUT_CC_MAG,              true,            Useless_bool,  Useless_bool   );
#  endif
#  ifdef GRAVITY
   ReadPara->Add( "OPT__OUTPUT_POT",            &OPT__OUTPUT_POT,                 false,           Useless_bool,  Useless_bool   );
#  endif
#  ifdef PARTICLE
   ReadPara->Add( "OPT__OUTPUT_PAR_DENS",       &OPT__OUTPUT_PAR_DENS,            PAR_OUTPUT_DENS_PAR_ONLY, 0,    2              );
#  endif
#  if ( MODEL == HYDRO )
   ReadPara->Add( "OPT__OUTPUT_PRES",           &OPT__OUTPUT_PRES,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_TEMP",           &OPT__OUTPUT_TEMP,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_CS",             &OPT__OUTPUT_CS,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_DIVVEL",         &OPT__OUTPUT_DIVVEL,              false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_MACH",           &OPT__OUTPUT_MACH,                false,           Useless_bool,  Useless_bool   );
#  endif
#  ifdef MHD
   ReadPara->Add( "OPT__OUTPUT_DIVMAG",         &OPT__OUTPUT_DIVMAG,              false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__OUTPUT_USER_FIELD",     &OPT__OUTPUT_USER_FIELD,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_MODE",           &OPT__OUTPUT_MODE,               -1,               1,             3              );
// do not check OUTPUT_STEP and OUTPUT_DT since they depend on OPT__OUTPUT_MODE
   ReadPara->Add( "OUTPUT_STEP",                &OUTPUT_STEP,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OUTPUT_DT",                  &OUTPUT_DT,                      -1.0,             NoMin_double,  NoMax_double   );
// do not check OUTPUT_PART_X/Y/Z since they depend on OPT__OUTPUT_PART
   ReadPara->Add( "OUTPUT_PART_X",              &OUTPUT_PART_X,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "OUTPUT_PART_Y",              &OUTPUT_PART_Y,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "OUTPUT_PART_Z",              &OUTPUT_PART_Z,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "INIT_DUMPID",                &INIT_DUMPID,                    -1,               NoMin_int,     NoMax_int      );


// yt inline analysis
#  ifdef SUPPORT_LIBYT
   ReadPara->Add( "YT_SCRIPT",                   YT_SCRIPT,                       Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "YT_VERBOSE",           (int*)&YT_VERBOSE,                      1,               0,             3              );
#  endif


// miscellaneous
   ReadPara->Add( "OPT__VERBOSE",               &OPT__VERBOSE,                    false,           Useless_bool,  Useless_bool   );
// do not check OPT__TIMING_BARRIER since it depends on other options
   ReadPara->Add( "OPT__TIMING_BARRIER",        &OPT__TIMING_BARRIER,            -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__TIMING_BALANCE",        &OPT__TIMING_BALANCE,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__TIMING_MPI",            &OPT__TIMING_MPI,                 false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_NOTE",           &OPT__RECORD_NOTE,                true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_UNPHY",          &OPT__RECORD_UNPHY,               true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_MEMORY",         &OPT__RECORD_MEMORY,              true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_PERFORMANCE",    &OPT__RECORD_PERFORMANCE,         true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__MANUAL_CONTROL",        &OPT__MANUAL_CONTROL,             true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_USER",           &OPT__RECORD_USER,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OPTIMIZE_AGGRESSIVE",   &OPT__OPTIMIZE_AGGRESSIVE,        false,           Useless_bool,  Useless_bool   );


// simulation checks
   ReadPara->Add( "OPT__CK_REFINE",             &OPT__CK_REFINE,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_PROPER_NESTING",     &OPT__CK_PROPER_NESTING,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_CONSERVATION",       &OPT__CK_CONSERVATION,            false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_NORMALIZE_PASSIVE",  &OPT__CK_NORMALIZE_PASSIVE,       false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_RESTRICT",           &OPT__CK_RESTRICT,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_FINITE",             &OPT__CK_FINITE,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_PATCH_ALLOCATE",     &OPT__CK_PATCH_ALLOCATE,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_FLUX_ALLOCATE",      &OPT__CK_FLUX_ALLOCATE,           false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO )
   ReadPara->Add( "OPT__CK_NEGATIVE",           &OPT__CK_NEGATIVE,                0,               0,             3              );
#  endif
   ReadPara->Add( "OPT__CK_MEMFREE",            &OPT__CK_MEMFREE,                 1.0,             0.0,           NoMax_double   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__CK_PARTICLE",           &OPT__CK_PARTICLE,                false,           Useless_bool,  Useless_bool   );
#  endif
#  ifdef MHD
   ReadPara->Add( "OPT__CK_INTERFACE_B",        &OPT__CK_INTERFACE_B,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_DIVERGENCE_B",       &OPT__CK_DIVERGENCE_B,            0,               0,             2              );
#  endif



// load parameters
// ********************************************************************************************************************************
   ReadPara->Read( FileName );
// ********************************************************************************************************************************


// free memory
   delete ReadPara;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Load_Parameter

