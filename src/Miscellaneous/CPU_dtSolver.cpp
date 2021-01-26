#include "GAMER.h"

#ifndef GPU


#if   ( MODEL == HYDRO )
void CPU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                            const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const int NPG,
                            const real dh, const real Safety, const real MinPres,
                            const EoS_DE2P_t EoS_DensEint2Pres_Func, const EoS_DP2C_t EoS_DensPres2CSqr_Func,
                            const EoS_TEM2C_t EoS_Temper2CSqr_Func, const EoS_GUESS_t EoS_GuessHTilde_Func,
                            const EoS_H2TEM_t EoS_HTilde2Temp_Func,
                            const double c_EoS_AuxArray_Flt[], const int c_EoS_AuxArray_Int[],
                            const real* const c_EoS_Table[EOS_NTABLE_MAX] );
#ifdef GRAVITY
void CPU_dtSolver_HydroGravity( real g_dt_Array[],
                                const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                const double g_Corner_Array[][3],
                                const int NPatchGroup, const real dh, const real Safety, const bool P5_Gradient,
                                const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                                const double c_ExtAcc_AuxArray[],
                                const double ExtAcc_Time );
#endif

#elif ( MODEL == ELBDM )

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver
// Description :  Invoke various CPU dt solvers
//
// Note        :  1. Invoked by InvokeSolver()
//                2. Corresponding preparation and closing steps are defined in dt_Prepare_XXX() and dt_Close()
//
// Parameter   :  TSolver      : Target dt solver
//                               --> DT_FLU_SOLVER : dt solver for fluid
//                                   DT_GRA_SOLVER : dt solver for gravity
//                dt_Array     : Array to store the minimum dt in each target patch
//                Flu_Array    : Array storing the prepared fluid     data of each target patch
//                Mag_Array    : Array storing the prepared B field   data of each target patch
//                Pot_Array    : Array storing the prepared potential data of each target patch
//                Corner_Array : Array storing the physical corner coordinates of each patch
//                NPatchGroup  : Number of patch groups evaluated simultaneously by GPU
//                dh           : Cell size
//                Safety       : dt safety factor
//                MinPres      : Minimum allowed pressure
//                P5_Gradient  : Use 5-points stencil to evaluate the potential gradient
//                UsePot       : Add self-gravity and/or external potential
//                ExtAcc       : Add external acceleration
//                TargetTime   : Target physical time
//
// Return      :  dt_Array
//-------------------------------------------------------------------------------------------------------
void CPU_dtSolver( const Solver_t TSolver, real dt_Array[], const real Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                   const real Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const real Pot_Array[][ CUBE(GRA_NXT) ],
                   const double Corner_Array[][3], const int NPatchGroup, const real dh, const real Safety,
                   const real MinPres, const bool P5_Gradient,
                   const bool UsePot, const OptExtAcc_t ExtAcc, const double TargetTime )
{

   switch ( TSolver )
   {
#     if   ( MODEL == HYDRO )
      case DT_FLU_SOLVER:
         CPU_dtSolver_HydroCFL( dt_Array, Flu_Array, Mag_Array, NPatchGroup, dh, Safety, MinPres,
                                EoS_DensEint2Pres_CPUPtr, EoS_DensPres2CSqr_CPUPtr, EoS_Temper2CSqr_CPUPtr,
                                EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      break;

#     ifdef GRAVITY
      case DT_GRA_SOLVER:
         CPU_dtSolver_HydroGravity( dt_Array, Pot_Array, Corner_Array, NPatchGroup, dh, Safety, P5_Gradient,
                                    UsePot, ExtAcc, CPUExtAcc_Ptr, ExtAcc_AuxArray, TargetTime );
      break;
#     endif


#     elif ( MODEL == ELBDM )


#     else
#        error : ERROR : unsupported MODEL !!
#     endif // MODEL

      default :
         Aux_Error( ERROR_INFO, "unsupported \"TSolver\" (%d) !!\n", TSolver );
   } // switch ( TSolver )

} // FUNCTION : CPU_dtSolver



#endif // #ifndef GPU
