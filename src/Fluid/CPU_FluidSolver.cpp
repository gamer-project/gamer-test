#ifndef GPU



#include "GAMER.h"
#include "CUFLU.h"

#ifndef GRAVITY
static double *ExtAcc_AuxArray = NULL;
static double *ExtPot_AuxArray = NULL;
static ExtAcc_t CPUExtAcc_Ptr  = NULL;
static ExtPot_t CPUExtPot_Ptr  = NULL;
#endif

#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
void CPU_FluidSolver_RTVD(
   real Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double Corner_Array[][3],
   const real Pot_Array_USG[][ CUBE(USG_NXT_F) ],
   const int NPatchGroup, const real dt, const real dh,
   const bool StoreFlux, const bool XYZ,
   const real MinDens, const real MinPres, const real MinEint,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func,
   const double c_EoS_AuxArray_Flt[],
   const int    c_EoS_AuxArray_Int[],
   const real* const c_EoS_Table[EOS_NTABLE_MAX] );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
void CPU_FluidSolver_MHM(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const int NPatchGroup,
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm,
   const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func,
   const EoS_GUESS_t EoS_GuessHTilde_Func,
   const EoS_H2TEM_t EoS_HTilde2Temp_Func,
   const EoS_TEM2H_t EoS_Temp2HTilde_Func,
   const EoS_TEM2C_t EoS_Temper2CSqr_Func,
   const double c_EoS_AuxArray_Flt[],
   const int    c_EoS_AuxArray_Int[],
   const real* const c_EoS_Table[EOS_NTABLE_MAX] );
#elif ( FLU_SCHEME == CTU )
void CPU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const int NPatchGroup, const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm,
   const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func,
   const double c_EoS_AuxArray_Flt[],
   const int    c_EoS_AuxArray_Int[],
   const real* const c_EoS_Table[EOS_NTABLE_MAX] );
#endif // FLU_SCHEME

#elif ( MODEL == ELBDM )
void CPU_ELBDMSolver( real Flu_Array_In [][FLU_NIN    ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real Flu_Array_Out[][FLU_NOUT   ][ PS2*PS2*PS2 ],
                      real Flux_Array[][9][NFLUX_TOTAL][ PS2*PS2 ],
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                      const real Taylor3_Coeff, const bool XYZ, const real MinDens );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL


#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*h_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*h_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*h_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*h_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*h_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*h_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#else
static real (*h_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
static real (*h_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ] = NULL;
#endif
#endif // FLU_SCHEME




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver
// Description :  1. MODEL == HYDRO : use CPU to solve the Euler equations by different schemes
//                2. MODEL == ELBDM : use CPU to solve the kinematic operator in the Schrodinger's equation
//
// Note        :  Currently five hydro schemes are supported in HYDRO :
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   3. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   4. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
//
// Parameter   :  h_Flu_Array_In      : Host array storing the input fluid variables
//                h_Flu_Array_Out     : Host array to store the output fluid variables
//                h_Mag_Array_In      : Host array storing the input B field (for MHD only)
//                h_Mag_Array_Out     : Host array to store the output B field (for MHD only)
//                h_DE_Array_Out      : Host array to store the dual-energy status
//                h_Flux_Array        : Host array to store the output fluxes (useful only if StoreFlux == true)
//                h_Ele_Array         : Host array to store the output electric field (for MHD only)
//                h_Corner_Array      : Host array storing the physical corner coordinates of each patch group
//                h_Pot_Array_USG     : Host array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup         : Number of patch groups to be evaluated
//                dt                  : Time interval to advance solution
//                dh                  : Grid size
//                StoreFlux           : true --> store the coarse-fine fluxes
//                StoreElectric       : true --> store the coarse-fine electric field
//                XYZ                 : true   : x->y->z ( forward sweep)
//                                      false1 : z->y->x (backward sweep)
//                                      --> only useful for the RTVD scheme
//                LR_Limiter          : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                      (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                     vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff        : Coefficient of the generalized MinMod limiter
//                ELBDM_Eta           : Particle mass / Planck constant
//                ELBDM_Taylor3_Coeff : Coefficient in front of the third term in the Taylor expansion for ELBDM
//                ELBDM_Taylor3_Auto  : true --> Determine ELBDM_Taylor3_Coeff automatically by invoking the
//                                               function "ELBDM_SetTaylor3Coeff"
//                Time                : Current physical time                      (for UNSPLIT_GRAVITY only)
//                UsePot              : Add self-gravity and/or external potential (for UNSPLIT_GRAVITY only)
//                ExtAcc              : Add external acceleration                  (for UNSPLIT_GRAVITY only)
//                MinDens/Pres/Eint   : Density, pressure, and internal energy floors
//                DualEnergySwitch    : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive         : true --> normalize passive scalars so that the sum of their mass density
//                                               is equal to the gas mass density
//                NNorm               : Number of passive scalars to be normalized
//                                      --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx             : Target variable indices to be normalized
//                                      --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres        : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff  : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver( real h_Flu_Array_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                      real h_Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                      real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                      real h_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                      char h_DE_Array_Out[][ CUBE(PS2) ],
                      real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                      real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                      const double h_Corner_Array[][3],
                      const real h_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
                      const int NPatchGroup, const real dt, const real dh,
                      const bool StoreFlux, const bool StoreElectric,
                      const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                      const real ELBDM_Eta, real ELBDM_Taylor3_Coeff, const bool ELBDM_Taylor3_Auto,
                      const double Time, const bool UsePot, const OptExtAcc_t ExtAcc,
                      const real MinDens, const real MinPres, const real MinEint,
                      const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                      const bool JeansMinPres, const real JeansMinPres_Coeff )
{

// check
#  ifdef GAMER_DEBUG
   if ( StoreFlux  &&  h_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "Flux_Array is not allocated !!\n" );

#  ifdef MHD
   if ( StoreElectric  &&  h_Ele_Array == NULL )   Aux_Error( ERROR_INFO, "Ele_Array is not allocated !!\n" );
#  endif

#  ifdef UNSPLIT_GRAVITY
   if ( UsePot  &&  h_Pot_Array_USG == NULL )   Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );

   if ( ExtAcc  &&  h_Corner_Array == NULL )    Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
#  endif
#  endif


#  if   ( MODEL == HYDRO )

#     if   ( FLU_SCHEME == RTVD )

      CPU_FluidSolver_RTVD( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, h_Corner_Array, h_Pot_Array_USG,
                            NPatchGroup, dt, dh, StoreFlux, XYZ, MinDens, MinPres, MinEint,
                            EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

#     elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

      CPU_FluidSolver_MHM ( h_Flu_Array_In, h_Flu_Array_Out, h_Mag_Array_In, h_Mag_Array_Out,
                            h_DE_Array_Out, h_Flux_Array, h_Ele_Array, h_Corner_Array, h_Pot_Array_USG,
                            h_PriVar, h_Slope_PPM, h_FC_Var, h_FC_Flux, h_FC_Mag_Half, h_EC_Ele,
                            NPatchGroup, dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff, Time,
                            UsePot, ExtAcc, CPUExtAcc_Ptr, ExtAcc_AuxArray, MinDens, MinPres, MinEint,
                            DualEnergySwitch, NormPassive, NNorm, NormIdx, JeansMinPres, JeansMinPres_Coeff,
                            EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                            EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_Temp2HTilde_CPUPtr,
                            EoS_Temper2CSqr_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

#     elif ( FLU_SCHEME == CTU )

      CPU_FluidSolver_CTU ( h_Flu_Array_In, h_Flu_Array_Out, h_Mag_Array_In, h_Mag_Array_Out,
                            h_DE_Array_Out, h_Flux_Array, h_Ele_Array, h_Corner_Array, h_Pot_Array_USG,
                            h_PriVar, h_Slope_PPM, h_FC_Var, h_FC_Flux, h_FC_Mag_Half, h_EC_Ele,
                            NPatchGroup, dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff, Time,
                            UsePot, ExtAcc, CPUExtAcc_Ptr, ExtAcc_AuxArray, MinDens, MinPres, MinEint,
                            DualEnergySwitch, NormPassive, NNorm, NormIdx, JeansMinPres, JeansMinPres_Coeff,
                            EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr,
                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

#     else

#     error : unsupported CPU hydro scheme

#     endif


#  elif ( MODEL == ELBDM )
//    evaluate the optimized Taylor expansion coefficient
      if ( ELBDM_Taylor3_Auto )  ELBDM_Taylor3_Coeff = ELBDM_SetTaylor3Coeff( dt, dh, ELBDM_Eta );

      CPU_ELBDMSolver( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, ELBDM_Eta, StoreFlux,
                       ELBDM_Taylor3_Coeff, XYZ, MinDens );

#  else
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL

} // FUNCTION : CPU_FluidSolver



#endif // #ifndef GPU
