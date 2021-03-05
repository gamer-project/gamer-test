#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Blast_Dens_Bg;       // background mass density
static double Blast_Pres_Bg;       // background pressure
static double Blast_Pres_Exp;      // explosion pressure
static double Blast_Radius;        // explosion radius
static double Blast_Center[3];     // explosion center
#ifdef MHD
static double Blast_BField;        // magnetic field strength along the diagonal direction
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

   if ( !OPT__INIT_RESTRICT )
      Aux_Error( ERROR_INFO, "OPT__INIT_RESTRICT must be enabled !!\n" );


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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Blast_Dens_Bg",     &Blast_Dens_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Bg",     &Blast_Pres_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Exp",    &Blast_Pres_Exp,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius",      &Blast_Radius,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Center_X",    &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   ReadPara->Add( "Blast_Center_Y",    &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   ReadPara->Add( "Blast_Center_Z",    &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );
#  ifdef MHD
   ReadPara->Add( "Blast_BField",      &Blast_BField,           5.0e-2,       NoMin_double,     NoMax_double      );
#  endif

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 5.0e-3;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
//    assuming EOS_GAMMA (must not invoke any EoS routine here since it has not been initialized)
      const double ExpVol  = 4.0*M_PI/3.0*CUBE(Blast_Radius);
      const double ExpEngy = Blast_Pres_Exp/(GAMMA-1.0)*ExpVol;
#     if ( EOS != EOS_GAMMA )
      Aux_Message( stderr, "WARNING : the total explosion energy below assumes EOS_GAMMA !!\n" );
#     endif

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background mass density   = %13.7e\n", Blast_Dens_Bg );
      Aux_Message( stdout, "  background pressure       = %13.7e\n", Blast_Pres_Bg );
      Aux_Message( stdout, "  explosion pressure        = %13.7e\n", Blast_Pres_Exp );
      Aux_Message( stdout, "  total explosion energy    = %13.7e (assuming constant-gamma EoS)\n", ExpEngy );
      Aux_Message( stdout, "  explosion radius          = %13.7e\n", Blast_Radius );
      Aux_Message( stdout, "  explosion center          = (%13.7e, %13.7e, %13.7e)\n", Blast_Center[0], Blast_Center[1],
                                                                                       Blast_Center[2] );
#     ifdef MHD
      Aux_Message( stdout, "  magnetic field strength   = %13.7e\n", Blast_BField );
#     endif
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

   const double r = SQRT( SQR(x-Blast_Center[0]) + SQR(y-Blast_Center[1]) + SQR(z-Blast_Center[2]) );
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = Blast_Dens_Bg;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Pres = ( r <= Blast_Radius ) ? Blast_Pres_Exp : Blast_Pres_Bg;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table, NULL ); // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );         // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

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

   magnetic[MAGX] = Blast_BField / sqrt(3.0);
   magnetic[MAGY] = Blast_BField / sqrt(3.0);
   magnetic[MAGZ] = Blast_BField / sqrt(3.0);

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_BlastWave()
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

} // FUNCTION : Init_TestProb_Hydro_BlastWave
