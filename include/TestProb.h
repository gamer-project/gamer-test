#ifndef __TESTPROB_H__
#define __TESTPROB_H__



// common headers
#include "ReadPara.h"


// common function prototypes
static void Validate();
static void SetParameter();
static void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                       const int lv, double AuxArray[] );
#ifdef MHD
static void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                         const int lv, double AuxArray[] );
#endif


// function pointers of various user-specified routines
extern void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] );
#ifdef MHD
extern void (*Init_Function_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                              const int lv, double AuxArray[] );
#endif
extern void (*Init_ByFile_User_Ptr)( real fluid_out[], const real fluid_in[], const int nvar_in,
                                     const double x, const double y, const double z, const double Time,
                                     const int lv, double AuxArray[] );
extern void (*Init_Field_User_Ptr)();
extern void (*Init_User_Ptr)();
extern void (*Output_User_Ptr)();
extern bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
extern void (*Mis_GetTimeStep_User_Ptr)( const int lv, const double dTime_dt );
extern void (*Aux_Record_User_Ptr)();
extern void (*BC_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                            const int lv, double AuxArray[] );
#ifdef MHD
extern void (*BC_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                   const int lv, double AuxArray[] );
#endif
extern bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] );
extern void (*End_User_Ptr)();
#ifdef GRAVITY
extern real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                               const int lv, double AuxArray[] );
extern void (*Poi_UserWorkBeforePoisson_Ptr)( const double Time, const int lv );
extern void (*Init_ExtAcc_Ptr)();
extern void (*Init_ExtPot_Ptr)();
extern void (*End_ExtAcc_Ptr)();
extern void (*End_ExtPot_Ptr)();
#endif
#ifdef PARTICLE
extern void (*Par_Init_ByFunction_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                        real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                        real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                        real *AllAttribute[PAR_NATT_TOTAL] );
extern void (*Par_Init_Attribute_User_Ptr)();
#endif
#if ( MODEL == HYDRO )
extern void (*EoS_Init_Ptr)();
extern void (*EoS_End_Ptr)();
#endif
extern void (*Src_Init_User_Ptr)();
extern void (*Init_DerivedField_User_Ptr)();
#ifdef FEEDBACK
extern void (*FB_Init_User_Ptr)();
#endif


// helper macro for printing warning messages when resetting parameters
#  define FORMAT_INT       %- 21d
#  define FORMAT_LONG      %- 21ld
#  define FORMAT_UINT      %- 21u
#  define FORMAT_ULONG     %- 21lu
#  define FORMAT_BOOL      %- 21d
#  define FORMAT_REAL      %- 21.14e
#  define PRINT_WARNING( name, var, format )                                                             \
   {                                                                                                     \
      if ( MPI_Rank == 0 )                                                                               \
         Aux_Message( stderr, "WARNING : parameter [%-25s] is reset to [" EXPAND_AND_QUOTE(format) "] "  \
                              "for the adopted test problem\n", name, var );                             \
   }



#endif // #ifndef __TESTPROB_H__
