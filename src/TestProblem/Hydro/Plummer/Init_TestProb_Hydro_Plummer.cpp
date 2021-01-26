#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
       int    Plummer_RSeed;        // random seed for setting particle position and velocity
       double Plummer_Rho0;         // peak density
       double Plummer_R0;           // scale radius
       double Plummer_MaxR;         // maximum radius for particles
       bool   Plummer_Collision;    // true/false --> two colliding Plummer clouds/single Plummer cloud
       double Plummer_Collision_D;  // distance between two colliding Plummer clouds
       double Plummer_Center[3];    // central coordinates
       double Plummer_BulkVel[3];   // bulk velocity
       double Plummer_GasMFrac;     // gas                   mass fraction
       double Plummer_ParMFrac;     // particle              mass fraction
       double Plummer_ExtAccMFrac;  // external acceleration mass fraction
       double Plummer_ExtPotMFrac;  // external potential    mass fraction
       int    Plummer_MassProfNBin; // number of radial bins in the mass profile table
static bool   Plummer_AddColor;     // assign different colors to different clouds for Plummer_Collision

static double Plummer_FreeT;        // free-fall time at Plummer_R0

static FieldIdx_t Plummer_Idx_Cloud0 = Idx_Undefined;    // field indices for Plummer_AddColor
static FieldIdx_t Plummer_Idx_Cloud1 = Idx_Undefined;
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_Plummer( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] );
#endif
void Init_ExtAcc_Plummer();
void Init_ExtPot_Plummer();




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
#  ifdef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be disabled !!\n" );
#  endif
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Plummer_RSeed",        &Plummer_RSeed,         123,           0,                NoMax_int         );
   ReadPara->Add( "Plummer_Rho0",         &Plummer_Rho0,          1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Plummer_R0",           &Plummer_R0,            0.1,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Plummer_MaxR",         &Plummer_MaxR,          0.375,         Eps_double,       NoMax_double      );
   ReadPara->Add( "Plummer_Collision",    &Plummer_Collision,     false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Plummer_Collision_D",  &Plummer_Collision_D,   1.5,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_CenterX",      &Plummer_Center[0],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_CenterY",      &Plummer_Center[1],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_CenterZ",      &Plummer_Center[2],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_BulkVelX",     &Plummer_BulkVel[0],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_BulkVelY",     &Plummer_BulkVel[1],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_BulkVelZ",     &Plummer_BulkVel[2],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Plummer_GasMFrac",     &Plummer_GasMFrac,      0.25,          Eps_double,       1.0               );
   ReadPara->Add( "Plummer_ExtAccMFrac",  &Plummer_ExtAccMFrac,   0.25,          0.0,              1.0               );
   ReadPara->Add( "Plummer_ExtPotMFrac",  &Plummer_ExtPotMFrac,   0.25,          0.0,              1.0               );
   ReadPara->Add( "Plummer_MassProfNBin", &Plummer_MassProfNBin,  1000,          2,                NoMax_int         );
   ReadPara->Add( "Plummer_AddColor",     &Plummer_AddColor,      false,         Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   for (int d=0; d<3; d++)
      if ( Plummer_Center[d] == NoDef_double )  Plummer_Center[d] = 0.5*amr->BoxSize[d];

   if ( !Plummer_Collision  &&  Plummer_AddColor )
      Aux_Error( ERROR_INFO, "\"Plummer_AddColor\" must work with \"Plummer_Collision\" !!\n" );

// (1-3) check and reset the runtime parameters
   if ( Plummer_AddColor  &&  NCOMP_PASSIVE_USER != 2 )
      Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_USER to 2 for \"Plummer_AddColor\" !!\n" );

#  ifdef GRAVITY
   if ( !OPT__EXT_ACC  &&  Plummer_ExtAccMFrac != 0.0 )
   {
      Plummer_ExtAccMFrac = 0.0;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : \"Plummer_ExtAccMFrac\" is reset to 0.0 since OPT__EXT_ACC is disabled !!\n" );
   }

   if ( !OPT__EXT_POT  &&  Plummer_ExtPotMFrac != 0.0 )
   {
      Plummer_ExtPotMFrac = 0.0;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : \"Plummer_ExtPotMFrac\" is reset to 0.0 since OPT__EXT_POT is disabled !!\n" );
   }
#  endif

#  ifndef PARTICLE
   Plummer_ParMFrac = 0.0;
#  endif

#  ifdef GRAVITY
   if ( OPT__SELF_GRAVITY )
   {
//    ensure the sum of all mass fractions equals unity
      const double NonParMFrac = Plummer_GasMFrac + Plummer_ExtAccMFrac + Plummer_ExtPotMFrac;

#     ifdef PARTICLE
      if ( NonParMFrac > 1.0 )
         Aux_Error( ERROR_INFO, "GasMFrac (%13.7e) + ExtAccMFrac (%13.7e) + ExtPotMFrac (%13.7e) = %13.7e > 1.0 !!\n",
                    Plummer_GasMFrac, Plummer_ExtAccMFrac, Plummer_ExtPotMFrac, NonParMFrac );

      else
         Plummer_ParMFrac = 1.0 - NonParMFrac;

#     else
	  if (  ! Mis_CompareRealValue( NonParMFrac, 1.0, NULL, false )  )
      {
         Plummer_GasMFrac = 1.0 - Plummer_ExtAccMFrac - Plummer_ExtPotMFrac;

         if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : \"Plummer_GasMFrac\" is reset to %13.7e !!\n", Plummer_GasMFrac );
      }
#     endif
   } // if ( OPT__SELF_GRAVITY )

   else
   {
//    without self-gravity, gas and particle mass normalization can be set arbitrarily as they contribute no gravity
//    --> just ensure the sum of their mass fractions equals unity
#     ifdef PARTICLE
      if ( Plummer_GasMFrac > 1.0 )
         Aux_Error( ERROR_INFO, "GasMFrac (%13.7e) > 1.0 !!\n", Plummer_GasMFrac );

      else
         Plummer_ParMFrac = 1.0 - Plummer_GasMFrac;

#     else
      if (  ! Mis_CompareRealValue( Plummer_GasMFrac, 1.0, NULL, false )  )
      {
         Plummer_GasMFrac = 1.0;

         if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : \"Plummer_GasMFrac\" is reset to %13.7e !!\n", Plummer_GasMFrac );
      }
#     endif

//    ensure the sum of the mass fractions of external acceleration and potential equals unity
      if ( OPT__EXT_POT )
      {
         if ( OPT__EXT_ACC )
         {
            if (  ! Mis_CompareRealValue( Plummer_ExtAccMFrac+Plummer_ExtPotMFrac, 1.0, NULL, false )  )
            {
               Plummer_ExtPotMFrac = 1.0 - Plummer_ExtAccMFrac;

               if ( MPI_Rank == 0 )
                  Aux_Message( stderr, "WARNING : \"Plummer_ExtPotMFrac\" is reset to %13.7e !!\n", Plummer_ExtPotMFrac );
            }
         }

         else
         {
            if (  ! Mis_CompareRealValue( Plummer_ExtPotMFrac, 1.0, NULL, false )  )
            {
               Plummer_ExtPotMFrac = 1.0;

               if ( MPI_Rank == 0 )
                  Aux_Message( stderr, "WARNING : \"Plummer_ExtPotMFrac\" is reset to %13.7e !!\n", Plummer_ExtPotMFrac );
            }
         }
      } // if ( OPT__EXT_POT )

      else
      {
         if ( OPT__EXT_ACC )
         {
            if (  ! Mis_CompareRealValue( Plummer_ExtAccMFrac, 1.0, NULL, false )  )
            {
               Plummer_ExtAccMFrac = 1.0;

               if ( MPI_Rank == 0 )
                  Aux_Message( stderr, "WARNING : \"Plummer_ExtAccMFrac\" is reset to %13.7e !!\n", Plummer_ExtAccMFrac );
            }
         }

         else
            Aux_Error( ERROR_INFO, "all gravity options are disabled (OPT__SELF_GRAVITY, OPT__EXT_ACC, OPT__EXT_POT) !!\n" );
      } // if ( OPT__EXT_POT ) ... else ...
   } // if ( OPT__SELF_GRAVITY ) ... else ...
#  endif // #ifdef GRAVITY

// (2) set the problem-specific derived parameters
#  ifdef GRAVITY
   Plummer_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Plummer_Rho0) );
#  endif


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = (Plummer_Collision) ? 50.0 : 20.0*Plummer_FreeT;

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
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     Plummer_RSeed );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", Plummer_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", Plummer_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", Plummer_MaxR );
      Aux_Message( stdout, "  test mode                                 = %s\n",    (Plummer_Collision)?
                                                                                     "colliding clouds":"single cloud" );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "  central coordinate [%d]                   = %14.7e\n", d, Plummer_Center[d] );
      if ( Plummer_Collision ) {
      Aux_Message( stdout, "  initial distance between two clouds       = %13.7e\n", Plummer_Collision_D );
      Aux_Message( stdout, "  assign colors to different clouds         = %d\n",     Plummer_AddColor ); }
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "  bulk velocity [%d]                        = %14.7e\n", d, Plummer_BulkVel[d] );
      Aux_Message( stdout, "  gas                   mass fraction       = %13.7e\n", Plummer_GasMFrac );
      Aux_Message( stdout, "  particle              mass fraction       = %13.7e\n", Plummer_ParMFrac );
      Aux_Message( stdout, "  external acceleration mass fraction       = %13.7e\n", Plummer_ExtAccMFrac );
      Aux_Message( stdout, "  external potential    mass fraction       = %13.7e\n", Plummer_ExtPotMFrac );
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     Plummer_MassProfNBin );
      Aux_Message( stdout, "  free-fall time at the scale radius        = %13.7e\n", Plummer_FreeT );
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

// gas share the same density profile as particles (except for different total masses)
   const double TotM    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
   const double GasRho0 = Plummer_Rho0*Plummer_GasMFrac;
   const double PresBg  = 0.0;   // background pressure (set to 0.0 by default)

#  ifdef SRHD
   real Prim[NCOMP_FLUID];
#  endif

   double Dens=0.0, MomX=0.0, MomY=0.0, MomZ=0.0, Pres=0.0, Eint, Etot;
   double r2, a2, Dens1Cloud;


   if ( Plummer_Collision )
   {
      const double Coll_Offset = 0.5*Plummer_Collision_D/sqrt(3.0);
      double Center[3];

      for (int t=-1; t<=1; t+=2)
      {
         for (int d=0; d<3; d++)    Center[d] = Plummer_Center[d] + Coll_Offset*(double)t;

         r2         = SQR(x-Center[0]) + SQR(y-Center[1]) + SQR(z-Center[2]);
         a2         = r2 / SQR(Plummer_R0);
         Dens1Cloud = GasRho0 * pow( 1.0 + a2, -2.5 );

         Dens += Dens1Cloud;
#        ifdef GRAVITY
         Pres += NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg;
#        endif

//       add different colors for different clouds
         if ( Plummer_AddColor )
         fluid[ (t==-1)?Plummer_Idx_Cloud0:Plummer_Idx_Cloud1 ] = Dens1Cloud;
      }

#     ifdef SRHD
      Prim[0] = Dens;
      Prim[1] = Plummer_BulkVel[0];
      Prim[2] = Plummer_BulkVel[1];
      Prim[3] = Plummer_BulkVel[2];
      Prim[4] = Pres;
#     else
      MomX = Dens*Plummer_BulkVel[0];
      MomY = Dens*Plummer_BulkVel[1];
      MomZ = Dens*Plummer_BulkVel[2];
#     endif
   } // if ( Plummer_Collision )

   else
   {
      r2         = SQR(x-Plummer_Center[0]) + SQR(y-Plummer_Center[1]) + SQR(z-Plummer_Center[2]);
      a2         = r2 / SQR(Plummer_R0);
      Dens1Cloud = GasRho0 * pow( 1.0 + a2, -2.5 );

#     ifdef GRAVITY
      Pres = NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg;
#     endif

#     ifdef SRHD
      Prim[0] = Dens1Cloud;
      Prim[1] = Plummer_BulkVel[0];
      Prim[2] = Plummer_BulkVel[1];
      Prim[3] = Plummer_BulkVel[2];
      Prim[4] = Pres;
#     else
      Dens = Dens1Cloud;
      MomX = Dens1Cloud*Plummer_BulkVel[0];
      MomY = Dens1Cloud*Plummer_BulkVel[1];
      MomZ = Dens1Cloud*Plummer_BulkVel[2];
#     endif


//    just set all passive scalars as zero
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;
   } // if ( Plummer_Collision ) ... else ...

#  ifdef SRHD
   Hydro_Pri2Con( Prim, fluid, NULL_BOOL, NULL_INT, NULL, NULL,
                  EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
#  else
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
#  endif


} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_Plummer
// Description :  Add the problem-specific fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_Plummer()
{

   if ( Plummer_AddColor )
   {
      Plummer_Idx_Cloud0 = AddField( "Cloud0", NORMALIZE_YES );
      Plummer_Idx_Cloud1 = AddField( "Cloud1", NORMALIZE_YES );
   }

} // FUNCTION : AddNewField_Plummer
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Plummer
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Plummer()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr  = SetGridIC;
   Init_Field_User_Ptr     = AddNewField_Plummer;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Plummer;
#  endif
#  ifdef GRAVITY
   Init_ExtAcc_Ptr         = Init_ExtAcc_Plummer;
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr         = Init_ExtPot_Plummer;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Plummer
