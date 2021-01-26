#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double ABC_Rho0;       // background density
static double ABC_P0;         // background pressure
static double ABC_V0;         // background velocity along the diagonal direction
static double ABC_CoeffA;     // coefficients A/B/C in the magnetic field IC (see Eq. [20] in Zhang et al., 2018, ApJS, 236, 50)
static double ABC_CoeffB;
static double ABC_CoeffC;
static int    ABC_NPeriod;    // number of periods along each direction
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


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifdef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be disabled !!!!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic !!\n" );


// warnings
   for (int s=0; s<6; s++)
      if ( OPT__BC_FLU[s] != BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : OPT__BC_FLU[%d] != BC_FLU_PERIODIC !?\n", s );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO && !defined SRHD )
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
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "ABC_Rho0",          &ABC_Rho0,              1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "ABC_P0",            &ABC_P0,                1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "ABC_V0",            &ABC_V0,                0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ABC_CoeffA",        &ABC_CoeffA,            1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ABC_CoeffB",        &ABC_CoeffB,            1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ABC_CoeffC",        &ABC_CoeffC,            1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ABC_NPeriod",       &ABC_NPeriod,           1,             1,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.5;

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
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background density  = %14.7e\n", ABC_Rho0    );
      Aux_Message( stdout, "  background pressure = %14.7e\n", ABC_P0      );
      Aux_Message( stdout, "  background velocity = %14.7e\n", ABC_V0      );
      Aux_Message( stdout, "  coefficient A       = %14.7e\n", ABC_CoeffA  );
      Aux_Message( stdout, "  coefficient B       = %14.7e\n", ABC_CoeffB  );
      Aux_Message( stdout, "  coefficient C       = %14.7e\n", ABC_CoeffC  );
      Aux_Message( stdout, "  number of periods   = %d\n",     ABC_NPeriod );
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
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = ABC_Rho0;
   MomX = ABC_Rho0*ABC_V0/sqrt(3.0);
   MomY = ABC_Rho0*ABC_V0/sqrt(3.0);
   MomZ = ABC_Rho0*ABC_V0/sqrt(3.0);
   Pres = ABC_P0;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

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

   const double k = 2.0*M_PI/amr->BoxSize[0]*ABC_NPeriod;   // assuming cubic domain

   magnetic[MAGX] = ABC_CoeffA*sin( k*z ) + ABC_CoeffC*cos( k*y );
   magnetic[MAGY] = ABC_CoeffB*sin( k*x ) + ABC_CoeffA*cos( k*z );
   magnetic[MAGZ] = ABC_CoeffC*sin( k*y ) + ABC_CoeffB*cos( k*x );

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_MHD_ABC
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_MHD_ABC()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO && !defined SRHD )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_MHD_ABC
