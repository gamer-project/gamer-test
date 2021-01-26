#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static int      Jet_NJet;                    // number of jets (1/2)
static double   Jet_BgDens;                  // ambient density
static double   Jet_BgTemp;                  // ambient temperature
static double   Jet_BgVel[3];                // ambient velocity
static double  *Jet_Radius        = NULL;    // radius of the cylinder-shape jet source
static double  *Jet_HalfHeight    = NULL;    // half height of the cylinder-shape jet source
static double  *Jet_SrcVel        = NULL;    // jet velocity
static double  *Jet_SrcDens       = NULL;    // jet density
static double  *Jet_SrcTemp       = NULL;    // jet temperature
static double (*Jet_Vec)[3]       = NULL;    // jet orientation vector (x,y,z) (NOT necessary to be a unit vector)
static double (*Jet_CenOffset)[3] = NULL;    // jet central coordinates offset

static double   Jet_BgPres;                  // ambient pressure
static double  *Jet_SrcPres       = NULL;    // jet pressure
static double (*Jet_Cen)[3]       = NULL;    // jet central coordinates
static double  *Jet_WaveK         = NULL;    // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double  *Jet_MaxDis        = NULL;    // maximum distance between the cylinder-shape jet source and the jet center
                                             // --> ( Jet_Radius^2 + Jet_HalfHeight^2 )^0.5
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

#  ifdef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be disabled !!!!\n" );
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
// load the number of jets first
   ReadPara->Add( "Jet_NJet",          &Jet_NJet,             -1,             1,                2                 );
   ReadPara->Read( FileName );
   delete ReadPara;
   ReadPara = new ReadPara_t; // for loading the remaining parameters

// allocate memory for the variables whose size depends on the number of jets
   Jet_Radius     = new double [Jet_NJet];
   Jet_HalfHeight = new double [Jet_NJet];
   Jet_SrcVel     = new double [Jet_NJet];
   Jet_SrcDens    = new double [Jet_NJet];
   Jet_SrcTemp    = new double [Jet_NJet];
   Jet_Vec        = new double [Jet_NJet][3];
   Jet_CenOffset  = new double [Jet_NJet][3];
   Jet_SrcPres    = new double [Jet_NJet];
   Jet_Cen        = new double [Jet_NJet][3];
   Jet_WaveK      = new double [Jet_NJet];
   Jet_MaxDis     = new double [Jet_NJet];

// load the input parameters independent of the number of jets
   ReadPara->Add( "Jet_BgDens",        &Jet_BgDens,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_BgTemp",        &Jet_BgTemp,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_BgVel_x",       &Jet_BgVel[0],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_y",       &Jet_BgVel[1],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_z",       &Jet_BgVel[2],          0.0,           NoMin_double,     NoMax_double      );

// load the input parameters of each jet
   ReadPara->Add( "Jet0_Radius",       &Jet_Radius    [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_HalfHeight",   &Jet_HalfHeight[0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcVel",       &Jet_SrcVel    [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcDens",      &Jet_SrcDens   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcTemp",      &Jet_SrcTemp   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_Vec_x",        &Jet_Vec       [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_y",        &Jet_Vec       [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_z",        &Jet_Vec       [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_x",  &Jet_CenOffset [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_y",  &Jet_CenOffset [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_z",  &Jet_CenOffset [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );

   if ( Jet_NJet == 2 ) {
   ReadPara->Add( "Jet1_Radius",       &Jet_Radius    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_HalfHeight",   &Jet_HalfHeight[1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcVel",       &Jet_SrcVel    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcDens",      &Jet_SrcDens   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcTemp",      &Jet_SrcTemp   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_Vec_x",        &Jet_Vec       [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_y",        &Jet_Vec       [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_z",        &Jet_Vec       [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_x",  &Jet_CenOffset [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_y",  &Jet_CenOffset [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_z",  &Jet_CenOffset [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   }

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) convert to code units
   Jet_BgDens           *= 1.0       / UNIT_D;
   Jet_BgTemp           *= Const_keV / UNIT_E;
   for (int d=0; d<3; d++)
   Jet_BgVel[d]         *= 1.0       / UNIT_V;

   for (int n=0; n<Jet_NJet; n++) {
   Jet_Radius    [n]    *= Const_kpc / UNIT_L;
   Jet_HalfHeight[n]    *= Const_kpc / UNIT_L;
   Jet_SrcVel    [n]    *= 1.0       / UNIT_V;
   Jet_SrcDens   [n]    *= 1.0       / UNIT_D;
   Jet_SrcTemp   [n]    *= Const_keV / UNIT_E;
   for (int d=0; d<3; d++)
   Jet_CenOffset [n][d] *= Const_kpc / UNIT_L;
   }


// (2) set the problem-specific derived parameters
   Jet_BgPres = Jet_BgDens*Jet_BgTemp/(MOLECULAR_WEIGHT*Const_amu/UNIT_M);

   for (int n=0; n<Jet_NJet; n++)
   {
      Jet_WaveK  [n] = 0.5*M_PI/Jet_HalfHeight[n];
      Jet_MaxDis [n] = sqrt( SQR(Jet_HalfHeight[n]) + SQR(Jet_Radius[n]) );
      Jet_SrcPres[n] = Jet_SrcDens[n]*Jet_SrcTemp[n]/(MOLECULAR_WEIGHT*Const_amu/UNIT_M);

      for (int d=0; d<3; d++)    Jet_Cen[n][d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[n][d];
   }


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
      Aux_Message( stdout, "  test problem ID     = %d\n",             TESTPROB_ID                           );
      Aux_Message( stdout, "  Jet_BgDens          = % 14.7e g/cm^3\n", Jet_BgDens*UNIT_D                     );
      Aux_Message( stdout, "  Jet_BgTemp          = % 14.7e keV\n",    Jet_BgTemp*UNIT_E/Const_keV           );
      Aux_Message( stdout, "  Jet_BgVel[x]        = % 14.7e cm/s\n",   Jet_BgVel[0]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_BgVel[y]        = % 14.7e cm/s\n",   Jet_BgVel[1]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_BgVel[z]        = % 14.7e cm/s\n",   Jet_BgVel[2]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_NJet            = %d\n",             Jet_NJet                              );
      for (int n=0; n<Jet_NJet; n++) {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Jet # %d:\n", n );
      Aux_Message( stdout, "     Jet_Radius       = % 14.7e kpc\n",    Jet_Radius    [n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "     Jet_HalfHeight   = % 14.7e kpc\n",    Jet_HalfHeight[n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "     Jet_SrcVel       = % 14.7e cm/s\n",   Jet_SrcVel    [n]*UNIT_V              );
      Aux_Message( stdout, "     Jet_SrcDens      = % 14.7e g/cm^3\n", Jet_SrcDens   [n]*UNIT_D              );
      Aux_Message( stdout, "     Jet_SrcTemp      = % 14.7e keV\n",    Jet_SrcTemp   [n]*UNIT_E/Const_keV    );
      Aux_Message( stdout, "     Jet_Vec[x]       = % 14.7e\n",        Jet_Vec       [n][0]                  );
      Aux_Message( stdout, "     Jet_Vec[y]       = % 14.7e\n",        Jet_Vec       [n][1]                  );
      Aux_Message( stdout, "     Jet_Vec[z]       = % 14.7e\n",        Jet_Vec       [n][2]                  );
      Aux_Message( stdout, "     Jet_CenOffset[x] = % 14.7e kpc\n",    Jet_CenOffset [n][0]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "     Jet_CenOffset[y] = % 14.7e kpc\n",    Jet_CenOffset [n][1]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "     Jet_CenOffset[z] = % 14.7e kpc\n",    Jet_CenOffset [n][2]*UNIT_L/Const_kpc ); }
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = Jet_BgDens;
   MomX = Dens*Jet_BgVel[0];
   MomY = Dens*Jet_BgVel[1];
   MomZ = Dens*Jet_BgVel[2];
   Pres = Jet_BgPres;
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



//-------------------------------------------------------------------------------------------------------
// Function    :  End_CollidingJets
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_CollidingJets()
{

   delete [] Jet_Radius;
   delete [] Jet_HalfHeight;
   delete [] Jet_SrcVel;
   delete [] Jet_SrcDens;
   delete [] Jet_SrcTemp;
   delete [] Jet_Vec;
   delete [] Jet_CenOffset;
   delete [] Jet_SrcPres;
   delete [] Jet_Cen;
   delete [] Jet_WaveK;
   delete [] Jet_MaxDis;

} // FUNCTION : End_CollidingJets



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_CollidingJets
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
bool Flu_ResetByUser_CollidingJets( real fluid[], const double x, const double y, const double z, const double Time,
                                    const int lv, double AuxArray[] )
{

   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   double TempVec[3], _Jet_VecAbs, Jet_SrcVel_xyz[3];
   real   MomSin;


// loop over all cells to add the jet source
   for (int n=0; n<Jet_NJet; n++)
   {
//    distance: jet center to mesh
      for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[n][d];
      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );


//    exclude cells far away from the jet source
      if ( Dis_c2m > Jet_MaxDis[n] )   continue;


//    vectors for calculating the distance between cells and the jet sources
      for (int d=0; d<3; d++)    TempVec[d] = Jet_Cen[n][d] + Jet_Vec[n][d];


//    velocity components along different directions
      _Jet_VecAbs = 1.0/sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );
      for (int d=0; d<3; d++)    Jet_SrcVel_xyz[d] = Jet_SrcVel[n] * Jet_Vec[n][d] * _Jet_VecAbs;


//    distance: temporary vector to mesh
      for (int d=0; d<3; d++)    Vec_v2m[d] = r[d] - TempVec[d];
      Dis_v2m = sqrt( SQR(Vec_v2m[0]) + SQR(Vec_v2m[1]) + SQR(Vec_v2m[2]) );


//    distance: jet center to temporary vector
      Dis_c2v = sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );


//    check whether or not the target cell is within the jet source
      S      = 0.5*( Dis_c2m + Dis_v2m + Dis_c2v );
      Area   = sqrt( S*(S-Dis_c2m)*(S-Dis_v2m)*(S-Dis_c2v) );
      Jet_dr = 2.0*Area/Dis_c2v;
      Jet_dh = sqrt( Dis_c2m*Dis_c2m - Jet_dr*Jet_dr );

      if ( Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
      {
//       reset the fluid variables within the jet source
         fluid[DENS] = Jet_SrcDens[n];

//       use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
         MomSin      = fluid[DENS]*sin( Jet_WaveK[n]*Jet_dh );
         MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[n][0] + Vec_c2m[1]*Jet_Vec[n][1] + Vec_c2m[2]*Jet_Vec[n][2] );
         fluid[MOMX] = MomSin*Jet_SrcVel_xyz[0];
         fluid[MOMY] = MomSin*Jet_SrcVel_xyz[1];
         fluid[MOMZ] = MomSin*Jet_SrcVel_xyz[2];

//       assuming EoS requires no passive scalars
         double Eint = EoS_DensPres2Eint_CPUPtr( Jet_SrcDens[n], Jet_SrcPres[n], NULL,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

//       do NOT include magnetic energy here
         fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint, 0.0 );

//       return immediately since we do NOT allow different jet source to overlap
         return true;

      } // if (  Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
   } // for (int n=0; n<Jet_NJet; n++)


   return false;

} // FUNCTION : Flu_ResetByUser_CollidingJets
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CollidingJets
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CollidingJets()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO && !defined SRHD )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_CollidingJets;
   End_User_Ptr             = End_CollidingJets;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CollidingJets
