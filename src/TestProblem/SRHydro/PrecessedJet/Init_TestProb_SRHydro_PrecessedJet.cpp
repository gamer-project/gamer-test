#include "GAMER.h"
#include "TestProb.h"


#if  ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================

static double   Jet_BgDens;                             // ambient density
static double   Jet_BgVel[3];                           // ambient 4-velocity
static double   Jet_Radius;                             // radius of the cylinder-shape jet source
static double   Jet_HalfHeight;                         // half height of the cylinder-shape jet source
static double   Jet_SrcVel;                             // jet 4-velocity
static double   Jet_SrcDens;                            // jet density
static double   Jet_Cone_Vec[3];                        // cone orientation vector (x,y,z) (NOT necessary to be a unit vector)
static double   Jet_CenOffset[3];                          // jet central coordinates offset

static double   Jet_BgTemp;                             // ambient temperature
static double   Jet_SrcTemp;                            // jet temperature
static double   Jet_Cen[3];                             // jet central coordinates
static double   Jet_WaveK;                              // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double   Jet_MaxDis;                             // maximum distance between the cylinder-shape jet source and the jet center
static double   Jet_Angular_Velocity;                   // precession angular velocity (degree per code_time)
static double   Jet_Angle;                              // precession angle in degree
static double   Jet_BurstStartTime;                     // start burst time in jet
static double   Jet_BurstEndTime;                       // end burst time in jet
static double   Jet_Burst4Vel;                          // burt 4-velocity
static double   Jet_BurstDens;                          // burst proper density
static double   Jet_BurstTemp;                          // burst temperature
static bool     Flag_Burst4Vel;
static bool     Flag_BurstDens;
static bool     Flag_BurstTemp;
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
#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL must be SR_HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



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

// declare the variables related to jet
   ReadPara->Add( "Jet_BgDens",             &Jet_BgDens,           -1.0,           Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_BgVel_x",            &Jet_BgVel[0],          0.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BgVel_y",            &Jet_BgVel[1],          0.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BgVel_z",            &Jet_BgVel[2],          0.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BgTemp",             &Jet_BgTemp,           -1.0,           Eps_double,     NoMax_double    );

   ReadPara->Add( "Jet_Radius",             &Jet_Radius    ,       -1.0,           Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfHeight",         &Jet_HalfHeight,       -1.0,           Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcVel",             &Jet_SrcVel    ,       -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_SrcDens",            &Jet_SrcDens   ,       -1.0,           Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcTemp",            &Jet_SrcTemp   ,       -1.0,           Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_Cone_Vec_x",         &Jet_Cone_Vec[0],       NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Cone_Vec_y",         &Jet_Cone_Vec[1],       NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Cone_Vec_z",         &Jet_Cone_Vec[2],       NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_x",        &Jet_CenOffset [0],     NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_y",        &Jet_CenOffset [1],     NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_z",        &Jet_CenOffset [2],     NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Angular_Velocity",   &Jet_Angular_Velocity,  NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Angle",              &Jet_Angle,             NoDef_double,  NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstStartTime",     &Jet_BurstStartTime,   -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstEndTime",       &Jet_BurstEndTime,     -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Burst4Vel",          &Jet_Burst4Vel,        -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstDens",          &Jet_BurstDens,        -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstTemp",          &Jet_BurstTemp,        -1.0,           NoMin_double,   NoMax_double    );
   ReadPara->Add( "Flag_Burst4Vel",         &Flag_Burst4Vel,        false,         Useless_bool,   Useless_bool    );
   ReadPara->Add( "Flag_BurstDens",         &Flag_BurstDens,        false,         Useless_bool,   Useless_bool    );
   ReadPara->Add( "Flag_BurstTemp",         &Flag_BurstTemp,        false,         Useless_bool,   Useless_bool    );
   

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-1) validation
   bool Burst = Flag_Burst4Vel | Flag_BurstDens | Flag_BurstTemp;

   if (  Burst &&  Jet_BurstEndTime <= Jet_BurstStartTime)
     Aux_Error( ERROR_INFO, "Jet_BurstEndTime <= Jet_BurstStartTime !!\n" );


   if ( Burst && Jet_BurstEndTime >= END_T )
     Aux_Error( ERROR_INFO, "Jet_BurstEndTime >= END_T !!\n" );


// (1-2) convert to code unit
   Jet_BgDens               *= 1.0       / UNIT_D;
   Jet_BgTemp               *= Const_keV / ( Const_mp * SQR(Const_c) );

   for (int d=0; d<3; d++)    
   Jet_BgVel[d]             *= Const_c   / UNIT_V;

   Jet_Radius               *= Const_kpc / UNIT_L;
   Jet_HalfHeight           *= Const_kpc / UNIT_L;
   Jet_SrcVel               *= Const_c   / UNIT_V;
   Jet_SrcTemp              *= Const_keV / ( Const_mp * SQR(Const_c) );
   Jet_SrcDens              *= 1.0       / UNIT_D;

   for (int d=0; d<3; d++)
   Jet_CenOffset[d]         *= Const_kpc / UNIT_L;

   Jet_Angular_Velocity     *= UNIT_T / ( 1e3*Const_yr );
   Jet_BurstStartTime       *= 1e3*Const_yr / UNIT_T;
   Jet_BurstEndTime         *= 1e3*Const_yr / UNIT_T;
   Jet_Burst4Vel            *= Const_c   / UNIT_V;
   Jet_BurstDens            *= 1.0       / UNIT_D;
   Jet_BurstTemp            *= Const_keV / ( Const_mp * SQR(Const_c) );


// (2) set the problem-specific derived parameters
      Jet_WaveK   = 0.5*M_PI/Jet_HalfHeight;
      Jet_MaxDis  = sqrt( SQR(Jet_HalfHeight) + SQR(Jet_Radius) );

      for (int d=0; d<3; d++)    Jet_Cen[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 50.0*Const_Myr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID      = %d\n",                 TESTPROB_ID                                   );
      Aux_Message( stdout, "  Jet_BgDens           = % 14.7e g/cm^3\n",     Jet_BgDens*UNIT_D                             );
      Aux_Message( stdout, "  Jet_BgTemp           = % 14.7e keV\n",        Jet_BgTemp*Const_mp*SQR(Const_c)/Const_keV    );
      Aux_Message( stdout, "  Jet_BgVel[x]         = % 14.7e c\n",          Jet_BgVel[0]                                  );
      Aux_Message( stdout, "  Jet_BgVel[y]         = % 14.7e c\n",          Jet_BgVel[1]                                  );
      Aux_Message( stdout, "  Jet_BgVel[z]         = % 14.7e c\n",          Jet_BgVel[2]                                  );
      Aux_Message( stdout, "  Jet_Radius           = % 14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc                   );
      Aux_Message( stdout, "  Jet_HalfHeight       = % 14.7e kpc\n",        Jet_HalfHeight*UNIT_L/Const_kpc               );
      Aux_Message( stdout, "  Jet_SrcVel           = % 14.7e c\n",          Jet_SrcVel                                    );
      Aux_Message( stdout, "  Jet_SrcDens          = % 14.7e g/cm^3\n",     Jet_SrcDens*UNIT_D                            );
      Aux_Message( stdout, "  Jet_SrcTemp          = % 14.7e keV\n",        Jet_SrcTemp*Const_mp*SQR(Const_c)/Const_keV   );
      Aux_Message( stdout, "  Jet_Cone_Vec[x]      = % 14.7e\n",            Jet_Cone_Vec[0]                               );
      Aux_Message( stdout, "  Jet_Cone_Vec[y]      = % 14.7e\n",            Jet_Cone_Vec[1]                               );
      Aux_Message( stdout, "  Jet_Cone_Vec[z]      = % 14.7e\n",            Jet_Cone_Vec[2]                               );
      Aux_Message( stdout, "  Jet_CenOffset[x]     = % 14.7e kpc\n",        Jet_CenOffset[0]*UNIT_L/Const_kpc             );
      Aux_Message( stdout, "  Jet_CenOffset[y]     = % 14.7e kpc\n",        Jet_CenOffset[1]*UNIT_L/Const_kpc             );
      Aux_Message( stdout, "  Jet_CenOffset[z]     = % 14.7e kpc\n",        Jet_CenOffset[2]*UNIT_L/Const_kpc             );
      Aux_Message( stdout, "  Jet_Angular_Velocity = % 14.7e degree/kyr\n", Jet_Angular_Velocity*1e3*Const_yr / UNIT_T    );
      Aux_Message( stdout, "  Jet_Angle            = % 14.7e\n",            Jet_Angle                                     );
      Aux_Message( stdout, "  Jet_BurstStartTime   = % 14.7e kyr \n",       Jet_BurstStartTime*UNIT_T/(1e3*Const_yr)      );
      Aux_Message( stdout, "  Jet_BurstEndTime     = % 14.7e kyr \n",       Jet_BurstEndTime*UNIT_T/(1e3*Const_yr)        );
      Aux_Message( stdout, "  Jet_Burst4Vel        = % 14.7e c \n",         Jet_Burst4Vel                                 );
      Aux_Message( stdout, "  Jet_BurstDens        = % 14.7e g/cm^3\n",     Jet_BurstDens*UNIT_D                          );
      Aux_Message( stdout, "  Jet_BurstTemp        = % 14.7e keV\n",        Jet_BurstTemp*Const_mp*SQR(Const_c)/Const_keV );
      Aux_Message( stdout, "  Flag_Burst4Vel       = % d\n",                Flag_Burst4Vel                                );
      Aux_Message( stdout, "  Flag_BurstDens       = % d\n",                Flag_BurstDens                                );
      Aux_Message( stdout, "  Flag_BurstTemp       = % d\n",                Flag_BurstTemp                                );
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
// 4-velocity
   double Pri4Vel[NCOMP_FLUID] = { Jet_BgDens, Jet_BgVel[0], Jet_BgVel[1], Jet_BgVel[2], Jet_BgTemp*Jet_BgDens };

// cast double to real
#  ifndef FLOAT8
   double Out[NCOMP_FLUID];

   SRHydro_Pri2Con(Pri4Vel, Out, GAMMA);

   fluid [0] = (real) Out[0];
   fluid [1] = (real) Out[1];
   fluid [2] = (real) Out[2];
   fluid [3] = (real) Out[3];
   fluid [4] = (real) Out[4];
#  else
   SRHydro_Pri2Con(Pri4Vel, fluid, GAMMA);
#  endif

} // FUNCTION : SetGridIC




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_PrecessedJet
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_PrecessedJet( real fluid[], const double x, const double y, const double z, const double Time,
                                    const int lv, double AuxArray[] )
{

   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   double Jet_SrcVel_xyz[3],Jet_BurstVel_xyz[3], Jet_Vec[3], Dis_Cone_Vec, Dis_Cone_Vec_2;
   double Cos_Phi, Sin_Phi, Sin_Angle, Cos_Angle, Sin_Omega, Cos_Omega;
   double Angle, Omega_t;
   double Pri4Vel[NCOMP_FLUID];
   real   MomSin;
   bool BurstDens_Yes, Burst4Vel_Yes, BurstTemp_Yes;

// distance: jet center to mesh
   for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[d];
   Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

// exclude cells far away from the jet source
   if ( Dis_c2m > Jet_MaxDis )   return false;

   Dis_Cone_Vec_2 = sqrt(SQR(Jet_Cone_Vec[0]) + SQR(Jet_Cone_Vec[1])); 
   Dis_Cone_Vec   = sqrt(SQR(Jet_Cone_Vec[0]) + SQR(Jet_Cone_Vec[1]) + SQR(Jet_Cone_Vec[2]));

   Angle   = Jet_Angle * M_PI / 180.0;
   Omega_t = Jet_Angular_Velocity * Time * M_PI / 180.0;

   Sin_Angle = sin(Angle);
   Cos_Angle = cos(Angle);

   Sin_Omega = sin(Omega_t);
   Cos_Omega = cos(Omega_t);
   
// coordinate transformation
   if ( Dis_Cone_Vec_2 >  __DBL_EPSILON__ )
   {
       Cos_Phi = Jet_Cone_Vec[0] / Dis_Cone_Vec_2;
       Sin_Phi = Jet_Cone_Vec[1] / Dis_Cone_Vec_2;

       Jet_Vec[0] = - Dis_Cone_Vec    * Sin_Phi * Sin_Angle * Cos_Omega
                    - Jet_Cone_Vec[2] * Cos_Phi * Sin_Angle * Sin_Omega
                    + Dis_Cone_Vec_2  * Cos_Phi * Cos_Angle;

       Jet_Vec[1] =   Dis_Cone_Vec    * Cos_Phi * Sin_Angle * Cos_Omega
                    - Jet_Cone_Vec[2] * Sin_Phi * Sin_Angle * Sin_Omega
                    + Dis_Cone_Vec_2  * Sin_Phi * Cos_Angle;

       Jet_Vec[2] =   Dis_Cone_Vec_2  * Sin_Angle * Sin_Omega 
                    + Jet_Cone_Vec[2] * Cos_Angle;

   }
   else
   {
      Jet_Vec[0] = Dis_Cone_Vec * Sin_Angle * Cos_Omega;
      Jet_Vec[1] = Dis_Cone_Vec * Sin_Angle * Sin_Omega;
      Jet_Vec[2] = Dis_Cone_Vec * Cos_Angle;
   }


//  distance: jet center to temporary vector
    Dis_c2v = sqrt( SQR(Jet_Vec[0]) + SQR(Jet_Vec[1]) + SQR(Jet_Vec[2]) );

//  velocity components along different directions
    for (int d=0; d<3; d++)
    {   
      Jet_SrcVel_xyz[d]   = Jet_SrcVel    * Jet_Vec[d] / Dis_c2v;
      Jet_BurstVel_xyz[d] = Jet_Burst4Vel * Jet_Vec[d] / Dis_c2v;
    }


//  distance: temporary vector to mesh
    for (int d=0; d<3; d++)    Vec_v2m[d] = r[d] - Jet_Cen[d] - Jet_Vec[d];
    Dis_v2m = sqrt( SQR(Vec_v2m[0]) + SQR(Vec_v2m[1]) + SQR(Vec_v2m[2]) );


//  check whether or not the target cell is within the jet source
    S      = 0.5*( Dis_c2m + Dis_v2m + Dis_c2v );
    Area   = sqrt( S*(S-Dis_c2m)*(S-Dis_v2m)*(S-Dis_c2v) );
    Jet_dr = 2.0*Area/Dis_c2v;
    Jet_dh = sqrt( Dis_c2m*Dis_c2m - Jet_dr*Jet_dr );

//  reset the fluid variables within the jet source
    if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )
    {
//     use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
       MomSin      = sin( Jet_WaveK*Jet_dh );
       MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[0] + Vec_c2m[1]*Jet_Vec[1] + Vec_c2m[2]*Jet_Vec[2] );

       BurstDens_Yes = ( Flag_BurstDens && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;
       Burst4Vel_Yes = ( Flag_Burst4Vel && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;
       BurstTemp_Yes = ( Flag_BurstTemp && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;

       Pri4Vel[0] = BurstDens_Yes ? Jet_BurstDens               : Jet_SrcDens;
       Pri4Vel[1] = Burst4Vel_Yes ? Jet_BurstVel_xyz[0]*MomSin  : Jet_SrcVel_xyz[0]*MomSin;
       Pri4Vel[2] = Burst4Vel_Yes ? Jet_BurstVel_xyz[1]*MomSin  : Jet_SrcVel_xyz[1]*MomSin;
       Pri4Vel[3] = Burst4Vel_Yes ? Jet_BurstVel_xyz[2]*MomSin  : Jet_SrcVel_xyz[2]*MomSin;
       Pri4Vel[4] = BurstTemp_Yes ? Jet_BurstTemp*Jet_BurstDens : Jet_SrcTemp*Jet_SrcDens;
        
//     cast double to real
#      ifndef FLOAT8
       double Out[NCOMP_FLUID];

       SRHydro_Pri2Con(Pri4Vel, Out, GAMMA);

       fluid [0] = (real) Out[0];
       fluid [1] = (real) Out[1];
       fluid [2] = (real) Out[2];
       fluid [3] = (real) Out[3];
       fluid [4] = (real) Out[4];
#      else
       SRHydro_Pri2Con(Pri4Vel, fluid, GAMMA);
#      endif

//     return immediately since we do NOT allow different jet source to overlap
       return true;

    } // if (  Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )


   return false;

} // FUNCTION : Flu_ResetByUser_Jets



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_PrecessedJet
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_PrecessedJet()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// validate the compilation flags and runtime parameters
   Validate();

// set the problem-specific runtime parameters
   SetParameter();

   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_PrecessedJet;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_Jets
#endif
