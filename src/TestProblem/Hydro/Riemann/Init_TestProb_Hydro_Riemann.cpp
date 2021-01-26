#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
typedef int Riemann_t;
const Riemann_t
   SOD_SHOCK_TUBE =  0
  ,STRONG_SHOCK   =  1
  ,TWO_SHOCKS     =  2
  ,EINFELDT_1203  =  3
  ,EINFELDT_1125  =  4
  ,SONIC_RARE     =  5
#ifdef MHD
  ,RJ2A           =  6
  ,TORRILHON      =  7
  ,BRIO_WU        =  8
#endif
  ,NOH            =  9
  ,USER_DEFINED   = 10
  ;

static Riemann_t Riemann_Prob;         // target Riemann problem
static int       Riemann_LR;           // wave propagation direction (>0/<0 --> positive/negative direction)
static int       Riemann_XYZ;          // wave propagation direction (0/1/2 --> x/y/z)
static double    Riemann_Pos;          // position of discontinuity
static double    Riemann_Width;        // width of discontinuity

static char      Riemann_Name[100];    // name of the target Riemann problem
static double    Riemann_RhoL;         // left-state density
static double    Riemann_VelL;         // left-state longitidual velocity
static double    Riemann_VelL_T1;      // left-state transverse velocity 1
static double    Riemann_VelL_T2;      // left-state transverse velocity 2
static double    Riemann_PreL;         // left-state pressure
static double    Riemann_RhoR;         // right-state density
static double    Riemann_VelR;         // right-state longitidual velocity
static double    Riemann_VelR_T1;      // right-state transverse velocity 1
static double    Riemann_VelR_T2;      // right-state transverse velocity 2
static double    Riemann_PreR;         // right-state pressure
static double    Riemann_EndT;         // end physical time
#ifdef MHD
static double    Riemann_Mag;          // longitidual B field
static double    Riemann_MagL_T1;      // left-state transverse B field 1
static double    Riemann_MagL_T2;      // left-state transverse B field 2
static double    Riemann_MagR_T1;      // right-state transverse B field 1
static double    Riemann_MagR_T2;      // right-state transverse B field 2
#endif
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Riemann_Prob",      &Riemann_Prob,          -1,            0,                11                );
   ReadPara->Add( "Riemann_LR",        &Riemann_LR,             1,            NoMin_int,        NoMax_int         );
   ReadPara->Add( "Riemann_XYZ",       &Riemann_XYZ,            0,            0,                2                 );
   ReadPara->Add( "Riemann_RhoL",      &Riemann_RhoL,           __DBL_MAX__,  __DBL_MIN__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_RhoR",      &Riemann_RhoR,           __DBL_MAX__,  __DBL_MIN__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelL",      &Riemann_VelL,           __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelR",      &Riemann_VelR,           __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_PreL",      &Riemann_PreL,           __DBL_MAX__,  __DBL_MIN__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_PreR",      &Riemann_PreR,           __DBL_MAX__,  __DBL_MIN__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelL_T1",   &Riemann_VelL_T1,        __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelL_T2",   &Riemann_VelL_T2,        __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelR_T1",   &Riemann_VelR_T1,        __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_VelR_T2",   &Riemann_VelR_T2,        __DBL_MAX__, -__DBL_MAX__,      __DBL_MAX__       );
   ReadPara->Add( "Riemann_Pos",       &Riemann_Pos,            NoDef_double, NoMin_double,     NoMax_double      );
   ReadPara->Add( "Riemann_Width",     &Riemann_Width,          NoDef_double, Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( Riemann_Pos   == NoDef_double )   Riemann_Pos   = amr->BoxCenter[Riemann_XYZ];
   if ( Riemann_Width == NoDef_double )   Riemann_Width = 1.0e-10/UNIT_L;  // mimic a step function

   switch ( Riemann_Prob )
   {
      case SOD_SHOCK_TUBE : Riemann_RhoL = 1.0;    Riemann_VelL = 0.0;  Riemann_PreL = 1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;  Riemann_PreR = 0.1;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.1; ;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "Sod's shock tube" );
                            break;

      case STRONG_SHOCK   : Riemann_RhoL = 1250.0;  Riemann_VelL = 0.0;  Riemann_PreL = 500.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR =  125.0;  Riemann_VelR = 0.0;  Riemann_PreR =   5.0;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.4;     Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "strong shock" );
                            break;

      case TWO_SHOCKS     : Riemann_RhoL = 1.0;  Riemann_VelL = 3.0;  Riemann_PreL = 1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 2.0;  Riemann_VelR = 1.0;  Riemann_PreR = 1.0;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.1;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "two shocks" );
                            break;

      case EINFELDT_1203  : Riemann_RhoL = 1.0;  Riemann_VelL = -2.0;  Riemann_PreL = GAMMA-1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +2.0;  Riemann_PreR = GAMMA-1.0;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.1;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "Einfeldt's 1-2-0-3" );
                            if ( GAMMA < 1.0 )  Aux_Error( ERROR_INFO, "GAMMA (%13.7e) < 1.0 !!\n", GAMMA );
                            break;

      case EINFELDT_1125  : Riemann_RhoL = 1.0;  Riemann_VelL = -1.0;  Riemann_PreL = 2.5*(GAMMA-1.0);  Riemann_VelL_T1 = -2.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +1.0;  Riemann_PreR = 2.5*(GAMMA-1.0);  Riemann_VelR_T1 = +2.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.1;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "Einfeldt's 1-1-2-5" );
                            if ( GAMMA < 1.0 )  Aux_Error( ERROR_INFO, "GAMMA (%13.7e) < 1.0 !!\n", GAMMA );
                            break;

      case SONIC_RARE     : Riemann_RhoL = 1.0;    Riemann_VelL = 0.75;  Riemann_PreL = 1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;   Riemann_PreR = 0.1;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.1;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "sonic rarefaction wave" );
                            break;

#     ifdef MHD
      case RJ2A           : Riemann_RhoL = 1.08;  Riemann_VelL = 1.2;  Riemann_PreL = 0.95;  Riemann_VelL_T1 = 0.01;  Riemann_VelL_T2 = 0.5;
                            Riemann_RhoR = 1.0;   Riemann_VelR = 0.0;  Riemann_PreR = 1.0;   Riemann_VelR_T1 = 0.0;   Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.2;  Riemann_Pos = 0.5;
                            Riemann_MagL_T1 = 3.6/sqrt(4.0*M_PI);  Riemann_MagL_T2 = 2.0/sqrt(4.0*M_PI);
                            Riemann_MagR_T1 = 4.0/sqrt(4.0*M_PI);  Riemann_MagR_T2 = 2.0/sqrt(4.0*M_PI);
                            Riemann_Mag     = 2.0/sqrt(4.0*M_PI);
                            sprintf( Riemann_Name, "RJ2a" );
                            break;

      case TORRILHON      : Riemann_RhoL = 1.0;  Riemann_VelL = 0.0;   Riemann_PreL = 1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 0.2;  Riemann_VelR = 0.0;   Riemann_PreR = 0.2;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.08;  Riemann_Pos = 0.5;
                            Riemann_MagL_T1 = 1.0;       Riemann_MagL_T2 = 0.0;
                            Riemann_MagR_T1 = cos(3.0);  Riemann_MagR_T2 = sin(3.0);
                            Riemann_Mag     = 1.0;
                            sprintf( Riemann_Name, "Torrilhon" );
                            break;

      case BRIO_WU        : Riemann_RhoL = 1.0;    Riemann_VelL = 0.0;   Riemann_PreL = 1.0;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;   Riemann_PreR = 0.1;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.08;  Riemann_Pos = 0.5;
                            Riemann_MagL_T1 = +1.0;  Riemann_MagL_T2 = 0.0;
                            Riemann_MagR_T1 = -1.0;  Riemann_MagR_T2 = 0.0;
                            Riemann_Mag     = 0.75;
                            sprintf( Riemann_Name, "Brio & Wu shock tube" );
                            break;
#     endif // #ifdef MHD

      case NOH            : Riemann_RhoL = 1.0;  Riemann_VelL = +1.0;  Riemann_PreL = 1.0e-6;  Riemann_VelL_T1 = 0.0;  Riemann_VelL_T2 = 0.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = -1.0;  Riemann_PreR = 1.0e-6;  Riemann_VelR_T1 = 0.0;  Riemann_VelR_T2 = 0.0;
                            Riemann_EndT = 0.5;  Riemann_Pos = 0.5;
#                           ifdef MHD
                            Riemann_Mag = Riemann_MagL_T1 = Riemann_MagL_T2 = Riemann_MagR_T1 = Riemann_MagR_T2 = 0.0;
#                           endif
                            sprintf( Riemann_Name, "Noh's strong shock" );
                            break;

      case USER_DEFINED   : sprintf( Riemann_Name, "user defined" );
                            break;

      default : Aux_Error( ERROR_INFO, "unsupported Riemann problem (%d) !!\n", Riemann_Prob );
   } // switch ( Riemann_Prob )

// (1-3) check the runtime parameters
   if ( Riemann_LR == 0 )  Aux_Error( ERROR_INFO, "Riemann_LR must not be zero !!\n" );

#  ifdef MHD
   if (  (int)Riemann_Prob != RJ2A  &&  (int)Riemann_Prob != TORRILHON  &&  (int)Riemann_Prob != BRIO_WU  )
      Aux_Message( stderr, "WARNING : B field is zero in the %s Riemann problem (Riemann_Prob = %d) !!\n",
                   Riemann_Name, Riemann_Prob );
#  endif


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = Riemann_EndT;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if (  ( Riemann_XYZ == 0 && OPT__OUTPUT_PART != OUTPUT_X )  ||
         ( Riemann_XYZ == 1 && OPT__OUTPUT_PART != OUTPUT_Y )  ||
         ( Riemann_XYZ == 2 && OPT__OUTPUT_PART != OUTPUT_Z )    )
   {
      OPT__OUTPUT_PART = ( Riemann_XYZ == 0 ) ? OUTPUT_X : ( Riemann_XYZ == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      PRINT_WARNING( "OPT__OUTPUT_PART", OPT__OUTPUT_PART, FORMAT_INT );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                   = %d\n",     TESTPROB_ID     );
      Aux_Message( stdout, "  target Riemann problem            = %s\n",     Riemann_Name    );
      Aux_Message( stdout, "  position of discontinuity         = %14.7e\n", Riemann_Pos     );
      Aux_Message( stdout, "  width    of discontinuity         = %14.7e\n", Riemann_Width   );
      Aux_Message( stdout, "  left-state density                = %14.7e\n", Riemann_RhoL    );
      Aux_Message( stdout, "  left-state longitudinal velocity  = %14.7e\n", Riemann_VelL    );
      Aux_Message( stdout, "  left-state transverse velocity 1  = %14.7e\n", Riemann_VelL_T1 );
      Aux_Message( stdout, "  left-state transverse velocity 2  = %14.7e\n", Riemann_VelL_T2 );
      Aux_Message( stdout, "  left-state pressure               = %14.7e\n", Riemann_PreL    );
      Aux_Message( stdout, "  right-state density               = %14.7e\n", Riemann_RhoR    );
      Aux_Message( stdout, "  right-state longitudinal velocity = %14.7e\n", Riemann_VelR    );
      Aux_Message( stdout, "  right-state transverse velocity 1 = %14.7e\n", Riemann_VelR_T1 );
      Aux_Message( stdout, "  right-state transverse velocity 2 = %14.7e\n", Riemann_VelR_T2 );
      Aux_Message( stdout, "  right-state pressure              = %14.7e\n", Riemann_PreR    );
#     ifdef MHD
      Aux_Message( stdout, "  longitudinal B field              = %14.7e\n", Riemann_Mag     );
      Aux_Message( stdout, "  left-state transverse B field 1   = %14.7e\n", Riemann_MagL_T1 );
      Aux_Message( stdout, "  left-state transverse B field 2   = %14.7e\n", Riemann_MagL_T2 );
      Aux_Message( stdout, "  right-state transverse B field 1  = %14.7e\n", Riemann_MagR_T1 );
      Aux_Message( stdout, "  right-state transverse B field 2  = %14.7e\n", Riemann_MagR_T2 );
#     endif
      Aux_Message( stdout, "  propagation direction             = %s%s\n",   ( Riemann_LR > 0 ) ? "+" : "-",
                                                                             ( Riemann_XYZ == 0 ) ? "x" :
                                                                             ( Riemann_XYZ == 1 ) ? "y" : "z" );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   double r, Pres, Eint;
   int    MomIdx[3];
#  ifdef SRHD
   real Prim[NCOMP_FLUID], Cons[NCOMP_FLUID];
#  endif

   switch ( Riemann_XYZ )
   {
      case 0 : r=x;  MomIdx[0]=MOMX;  MomIdx[1]=MOMY;  MomIdx[2]=MOMZ;  break;
      case 1 : r=y;  MomIdx[0]=MOMY;  MomIdx[1]=MOMZ;  MomIdx[2]=MOMX;  break;
      case 2 : r=z;  MomIdx[0]=MOMZ;  MomIdx[1]=MOMX;  MomIdx[2]=MOMY;  break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Riemann_XYZ", Riemann_XYZ );
   }

   const double ds   = ( r - Riemann_Pos ) / Riemann_Width;
   const double Tanh = tanh( ds )*SIGN( Riemann_LR );
   const double dRho = 0.5*( Riemann_RhoR    - Riemann_RhoL    );
   const double aRho = 0.5*( Riemann_RhoR    + Riemann_RhoL    );
   const double dVel = 0.5*( Riemann_VelR    - Riemann_VelL    );
   const double aVel = 0.5*( Riemann_VelR    + Riemann_VelL    );
   const double dVT1 = 0.5*( Riemann_VelR_T1 - Riemann_VelL_T1 );
   const double aVT1 = 0.5*( Riemann_VelR_T1 + Riemann_VelL_T1 );
   const double dVT2 = 0.5*( Riemann_VelR_T2 - Riemann_VelL_T2 );
   const double aVT2 = 0.5*( Riemann_VelR_T2 + Riemann_VelL_T2 );
   const double dPre = 0.5*( Riemann_PreR    - Riemann_PreL    );
   const double aPre = 0.5*( Riemann_PreR    + Riemann_PreL    );

#  ifdef SRHD
   Prim[0] = (real)aRho + dRho*Tanh;
   Prim[1] = (real)aVel + dVel*Tanh; 
   Prim[2] = (real)aVT1 + dVT1*Tanh;
   Prim[3] = (real)aVT2 + dVT2*Tanh;
   Prim[4] = (real)aPre + dPre*Tanh;

   Hydro_Pri2Con( Prim, fluid, NULL_BOOL, NULL_INT, NULL, NULL,
                  EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                  EoS_AuxArray_Int, h_EoS_Table, NULL );
#  else
   fluid[ DENS      ] =   aRho + dRho*Tanh;
   fluid[ MomIdx[0] ] = ( aVel + dVel*Tanh )*fluid[DENS];
   fluid[ MomIdx[1] ] = ( aVT1 + dVT1*Tanh )*fluid[DENS];
   fluid[ MomIdx[2] ] = ( aVT2 + dVT2*Tanh )*fluid[DENS];
   Pres               =   aPre + dPre*Tanh;
#  endif

   if ( Riemann_LR < 0 )
   {
      fluid[ MomIdx[0] ] *= -1.0;
      fluid[ MomIdx[1] ] *= -1.0;
      fluid[ MomIdx[2] ] *= -1.0;
   }

#ifndef SRHD
// compute and store the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( fluid[DENS], Pres, fluid+NCOMP_FLUID,
                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

// do NOT include magnetic energy here
   fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint, 0.0 );
#endif

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   double r;
   int    DirL, DirT1, DirT2;

// determine the longitudinal and transverse directions
   switch ( Riemann_XYZ )
   {
      case 0 : r=x;  DirL=MAGX;  DirT1=MAGY;  DirT2=MAGZ;  break;
      case 1 : r=y;  DirL=MAGY;  DirT1=MAGZ;  DirT2=MAGX;  break;
      case 2 : r=z;  DirL=MAGZ;  DirT1=MAGX;  DirT2=MAGY;  break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Riemann_XYZ", Riemann_XYZ );
   }


// set B field
   const double ds      = ( r - Riemann_Pos ) / Riemann_Width;
   const double Tanh    = tanh( ds )*SIGN( Riemann_LR );
   const double dMag_T1 = 0.5*( Riemann_MagR_T1 - Riemann_MagL_T1 );
   const double aMag_T1 = 0.5*( Riemann_MagR_T1 + Riemann_MagL_T1 );
   const double dMag_T2 = 0.5*( Riemann_MagR_T2 - Riemann_MagL_T2 );
   const double aMag_T2 = 0.5*( Riemann_MagR_T2 + Riemann_MagL_T2 );

// longitudinal component
   magnetic[DirL ] = Riemann_Mag;

// transverse components
   magnetic[DirT1] = aMag_T1 + dMag_T1*Tanh;
   magnetic[DirT2] = aMag_T2 + dMag_T2*Tanh;


// change the B field sign if wave propagates along the negative direction
   if ( Riemann_LR < 0 )
      for (int v=0; v<NCOMP_MAG; v++)  magnetic[v] *= -1.0;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Riemann
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Riemann()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Riemann
