#include "GAMER.h"
#include "TestProb.h"

#include"Particle_IC_Constructor.h"
#include"string"//NEW INCLUDES

using namespace std;
extern Particle_IC_Constructor constructor_Models;//First Declared in Par_Init_ByFunction_Models
static double NFW_FreeT=10;

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_Models( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] );
#endif





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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetTotalTime_NFW
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
void SetTotalTime_NFW()
{
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  20.0*NFW_FreeT;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


} // FUNCTION : SetTotalTime_NFW
#endif


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
   double NFW_R0 = 0.1;
   double NFW_Rho0 = 1.0;
   double NFW_GasMFrac =0.001;
   double NFW_Center[3] = {1.5,1.5,1.5};

// gas share the same density profile as particles (except for different total masses)
   const double TotM    = 4.0/3.0*M_PI*CUBE(NFW_R0)*NFW_Rho0;
   const double GasRho0 = NFW_Rho0*NFW_GasMFrac;
   const double PresBg  = 0.0;   // background pressure (set to 0.0 by default)

   double r2, a2, Dens;


   
   {
      r2   = SQR(x-NFW_Center[0]) + SQR(y-NFW_Center[1]) + SQR(z-NFW_Center[2]);
      a2   = r2 / SQR(NFW_R0);
      Dens = GasRho0 * pow( 1.0 + a2, -2.5 );

      fluid[DENS] = Dens;
      fluid[MOMX] = 0;
      fluid[MOMY] = 0;
      fluid[MOMZ] = 0;
#     ifdef GRAVITY
      fluid[ENGY] = (  NEWTON_G*TotM*GasRho0 / ( 6.0*NFW_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
#     endif

//    just set all passive scalars as zero
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;
   } 


} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_NFW
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_NFW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetTotalTime_NFW();

   // User input
   bool single_cloud = false;
   //You may observe different scenarios if you choose single_cloud to be true or false.
   //If true, you may observe an NFW cloud evolving stably with time.
   //If false, you may observe two neighbouring NFW clouds torn apart by each other's tidal force. 
   
   if(single_cloud)
   {
      int NFW_num = 1;//User Input
      //The scenario is a single NFW cloud

      //User Input
      vector <string> TestProb_FileName;     //Test problem input parameter filenames
      vector <string> TypeName;                               //Type of models 
      vector <string> Profile_FileName;
      TestProb_FileName.push_back("Input__TestProb1");
      TypeName.push_back("NFW");
      Profile_FileName.push_back("NONE");

      //"NONE" means no additinal profile table.
      constructor_Models.construct_ic(NFW_num,TestProb_FileName,TypeName,Profile_FileName);
   }
   


   else
   {
      int NFW_num = 2;//User Input
      //The scenario is two NFW clouds, one is generated by NFW density function, the other is by a density profile table.

      //User Input
      vector <string> TestProb_FileName;     //Test problem input parameter filenames
      vector <string> TypeName;                               //Type of models 
      vector <string> Profile_FileName;
      TestProb_FileName.push_back("Input__TestProb1");
      TestProb_FileName.push_back("Input__TestProb2");
      TypeName.push_back("UNKNOWN");
      TypeName.push_back("Plummer");
      Profile_FileName.push_back("profile.txt");
      Profile_FileName.push_back("NONE");
      //"UNKNOWM" means generating paritcles using density profile table.
      //"NONE" means no additinal profile table.
      //"profile.txt" is the profile table file name.
      constructor_Models.construct_ic(NFW_num,TestProb_FileName,TypeName,Profile_FileName);
   }

   Init_Function_User_Ptr  = SetGridIC;
   
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Models;
#  endif

#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_NFW
