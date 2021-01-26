#include "GAMER.h"
#include "TestProb.h"

static double RandomNumber(RandomNumber_t *RNG, const double Min, const double Max );



// problem-specific global variables
// =======================================================================================
static int    KH_RSeed;             // random seed
static double KH_RAmp;              // random number amplitude in both vx, vy, and vz
static double KH_Pres;              // background pressure
static double KH_Vx1;               // background velocity x in the upper region
static double KH_Vx2;               // background velocity x in the lower region
static double KH_Vy1;               // background velocity y in the upper region
static double KH_Vy2;               // background velocity y in the lower region
static double KH_Rho1;              // mass density in the upper region
static double KH_Rho2;              // mass density in the lower region
       int    KH_RefineShearMaxLv;  // refine the regions around the shear plane to this level
       int    KH_PeriodicZFactor;   // decompose the simulation domain further into N periodic subdomains along the z axis

static RandomNumber_t *RNG = NULL;
// =======================================================================================

// problem-specific function prototypes
bool Flag_KelvinHelmholtzInstability( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );




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

#  ifdef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be disabled !!!!\n" );
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

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 1\" (i.e., periodic BC) !!\n" );

   if ( OPT__INIT_GRID_WITH_OMP )
      Aux_Error( ERROR_INFO, "please disable \"OPT__INIT_GRID_WITH_OMP\" !!\n" );


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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",     &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "KH_RSeed",            &KH_RSeed,              -1,            0,                NoMax_int         );
   ReadPara->Add( "KH_RAmp",             &KH_RAmp,               -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "KH_Pres",             &KH_Pres,               -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "KH_Vx1",              &KH_Vx1,                NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "KH_Vx2",              &KH_Vx2,                NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "KH_Vy1",              &KH_Vy1,                NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "KH_Vy2",              &KH_Vy2,                NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "KH_Rho1",             &KH_Rho1,               -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "KH_Rho2",             &KH_Rho2,               -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "KH_RefineShearMaxLv", &KH_RefineShearMaxLv,   2,             0,                TOP_LEVEL         );
   ReadPara->Add( "KH_PeriodicZFactor",  &KH_PeriodicZFactor,    1,             1,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// get the number of OpenMP threads
   int NT;
#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {  NT = omp_get_num_threads();  }
#  else
   {  NT = 1;                      }
#  endif

// allocate RNG
   RNG = new RandomNumber_t( NT );

// set the random seed of each MPI rank
   for (int t=0; t<NT; t++) {
      RNG->SetSeed(t, KH_RSeed+MPI_Rank*1000+t);
   }

// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 1.0;
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
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",      TESTPROB_ID );
      Aux_Message( stdout, "  KH_RSeed            = %d\n",      KH_RSeed );
      Aux_Message( stdout, "  KH_RAmp             = % 14.7e\n", KH_RAmp );
      Aux_Message( stdout, "  KH_Pres             = % 14.7e\n", KH_Pres );
      Aux_Message( stdout, "  KH_Vx1              = % 14.7e\n", KH_Vx1 );
      Aux_Message( stdout, "  KH_Vx2              = % 14.7e\n", KH_Vx2 );
      Aux_Message( stdout, "  KH_Vy1              = % 14.7e\n", KH_Vy1 );
      Aux_Message( stdout, "  KH_Vy2              = % 14.7e\n", KH_Vy2 );
      Aux_Message( stdout, "  KH_Rho1             = % 14.7e\n", KH_Rho1 );
      Aux_Message( stdout, "  KH_Rho2             = % 14.7e\n", KH_Rho2 );
      Aux_Message( stdout, "  KH_RefineShearMaxLv = %d\n",      KH_RefineShearMaxLv );
      Aux_Message( stdout, "  KH_PeriodicZFactor  = %d\n",      KH_PeriodicZFactor );
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

   const double dz_periodic = amr->BoxSize[2] / KH_PeriodicZFactor;
   const double z_shear     = 0.5*dz_periodic;
   const double z_periodic  = fmod( z, dz_periodic );

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

// region 1
   if ( z_periodic >= z_shear )
   {
      Dens = KH_Rho1;
      MomX = KH_Rho1*( KH_Vx1 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      MomY = KH_Rho1*( KH_Vy1 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      MomZ = KH_Rho1*(    0.0 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      Pres = KH_Pres;
   }

// region 2
   else
   {
      Dens = KH_Rho2;
      MomX = KH_Rho2*( KH_Vx2 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      MomY = KH_Rho2*( KH_Vy2 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      MomZ = KH_Rho2*(    0.0 + RandomNumber(RNG,-KH_RAmp,KH_RAmp) );
      Pres = KH_Pres;
   }

// compute the total gas energy
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
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  RandomNumber
// Description :  Generate a single random number
//
// Parameter   :  RNG : Random number generator
//                Min : Lower limit of the random number
//                Max : Upper limit of the random number
//
// Return      :  Random number
//-------------------------------------------------------------------------------------------------------
double RandomNumber(RandomNumber_t *RNG, const double Min, const double Max )
{

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   return RNG->GetValue( TID, Min, Max );

} // FUNCTION : RandomNumber



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_KelvinHelmholtzInstability
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_KelvinHelmholtzInstability()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO && !defined SRHD )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
   Flag_User_Ptr          = Flag_KelvinHelmholtzInstability;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_KelvinHelmholtzInstability
