#include "GAMER.h"
#include "TestProb.h"

#include <string>

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

// problem-specific global variables
// =======================================================================================
       int     Merger_Coll_NumHalos;      // number of clusters
static char    Merger_File_Prof1[1000];   // profile table of cluster 1
static char    Merger_File_Prof2[1000];   // profile table of cluster 2
static char    Merger_File_Prof3[1000];   // profile table of cluster 3
       char    Merger_File_Par1 [1000];   // particle file of cluster 1
       char    Merger_File_Par2 [1000];   // particle file of cluster 2
       char    Merger_File_Par3 [1000];   // particle file of cluster 3
       bool    Merger_Coll_IsGas1;        // (true/false) --> does cluster 1 have gas
       bool    Merger_Coll_IsGas2;        // (true/false) --> does cluster 2 have gas
       bool    Merger_Coll_IsGas3;        // (true/false) --> does cluster 3 have gas
       bool    Merger_Coll_UseMetals;     // (true/false) --> do the clusters have a metal field
       double  Merger_Coll_PosX1;         // x-position of the first cluster
       double  Merger_Coll_PosY1;         // y-position of the first cluster
       double  Merger_Coll_PosX2;         // x-position of the second cluster
       double  Merger_Coll_PosY2;         // y-position of the second cluster
       double  Merger_Coll_PosX3;         // x-position of the third cluster
       double  Merger_Coll_PosY3;         // y-position of the third cluster
       double  Merger_Coll_VelX1;         // x-velocity of the first cluster
       double  Merger_Coll_VelY1;         // y-velocity of the first cluster
       double  Merger_Coll_VelX2;         // x-velocity of the second cluster
       double  Merger_Coll_VelY2;         // y-velocity of the second cluster
       double  Merger_Coll_VelX3;         // x-velocity of the third cluster
       double  Merger_Coll_VelY3;         // y-velocity of the third cluster
       double  Merger_Coll_ColorRad1;     // "color" radius of the first cluster
       double  Merger_Coll_ColorRad2;     // "color" radius of the second cluster
       double  Merger_Coll_ColorRad3;     // "color" radius of the third cluster

static double *Table_R1 = NULL;           // radius of cluster 1
static double *Table_D1 = NULL;           // density of cluster 1
static double *Table_P1 = NULL;           // pressure of cluster 1
static double *Table_M1 = NULL;           // metallicity of cluster 1

static double *Table_R2 = NULL;           // radius of cluster 2
static double *Table_D2 = NULL;           // density of cluster 2
static double *Table_P2 = NULL;           // pressure of cluster 2
static double *Table_M2 = NULL;           // metallicity of cluster 2

static double *Table_R3 = NULL;           // radius of cluster 3
static double *Table_D3 = NULL;           // density of cluster 3
static double *Table_P3 = NULL;           // pressure of cluster 3
static double *Table_M3 = NULL;           // metallicity of cluster 3

static int     Merger_NBin1;              // number of radial bins of cluster 1
static int     Merger_NBin2;              // number of radial bins of cluster 2
static int     Merger_NBin3;              // number of radial bins of cluster 3

       FieldIdx_t ParTypeIdx     = Idx_Undefined;
static FieldIdx_t ColorField1Idx = Idx_Undefined;
static FieldIdx_t ColorField2Idx = Idx_Undefined;
static FieldIdx_t ColorField3Idx = Idx_Undefined;

// =======================================================================================

// problem-specific function prototypes
#ifdef MASSSIVE_PARTICLES
void Par_Init_ByFunction_ClusterMerger(const long NPar_ThisRank,
                                       const long NPar_AllRank,
                                       real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                       real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                       real *ParType, real *AllAttribute[PAR_NATT_TOTAL]);
#endif

int Read_Num_Points_ClusterMerger(std::string filename);
void Read_Profile_ClusterMerger(std::string filename, std::string fieldname,
                                double field[]);
void AddNewField_ClusterMerger();

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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef MASSIVE_PARTICLES
   Aux_Error( ERROR_INFO, "MASSIVE_PARTICLES must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_FLU* = 1) for this test !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_POT = 1) for this test !!\n" );
#  endif

#  ifdef MASSIVE_PARTICLES
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
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

#ifdef SUPPORT_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );

// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",        &VARIABLE,               DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************
   ReadPara->Add( "Merger_Coll_NumHalos",   &Merger_Coll_NumHalos,   2,                1,             3              );
   ReadPara->Add( "Merger_Coll_IsGas1",     &Merger_Coll_IsGas1,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_IsGas2",     &Merger_Coll_IsGas2,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_IsGas3",     &Merger_Coll_IsGas3,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_File_Prof1",       Merger_File_Prof1,      Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par1",        Merger_File_Par1,       Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Prof2",       Merger_File_Prof2,      Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par2",        Merger_File_Par2,       Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Prof3",       Merger_File_Prof3,      Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par3",        Merger_File_Par3,       Useless_str,      Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_Coll_PosX1",      &Merger_Coll_PosX1,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY1",      &Merger_Coll_PosY1,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosX2",      &Merger_Coll_PosX2,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY2",      &Merger_Coll_PosY2,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosX3",      &Merger_Coll_PosX3,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY3",      &Merger_Coll_PosY3,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX1",      &Merger_Coll_VelX1,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY1",      &Merger_Coll_VelY1,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX2",      &Merger_Coll_VelX2,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY2",      &Merger_Coll_VelY2,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX3",      &Merger_Coll_VelX3,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY3",      &Merger_Coll_VelY3,      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_UseMetals",  &Merger_Coll_UseMetals,  true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_ColorRad1",  &Merger_Coll_ColorRad1,  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_ColorRad2",  &Merger_Coll_ColorRad2,  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_ColorRad3",  &Merger_Coll_ColorRad3,  -1.0,             NoMin_double,  NoMax_double   );

   ReadPara->Read( FileName );

   delete ReadPara;

// Validate that we have the correct number of passive scalars

   if ( Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals != NCOMP_PASSIVE_USER )
      Aux_Error( ERROR_INFO,
                 "please set NCOMP_PASSIVE_USER (currently %d) == Merger_Coll_NumHalos + Merger_Coll_UseMetals (currently %d) in the Makefile !!\n", NCOMP_PASSIVE_USER, Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals );

// convert to code units
   Merger_Coll_PosX1 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY1 *= Const_kpc / UNIT_L;
   Merger_Coll_PosX2 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY2 *= Const_kpc / UNIT_L;
   Merger_Coll_PosX3 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY3 *= Const_kpc / UNIT_L;
   Merger_Coll_VelX1 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY1 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelX2 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY2 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelX3 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY3 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_ColorRad1 *= Const_kpc / UNIT_L;
   Merger_Coll_ColorRad2 *= Const_kpc / UNIT_L;
   Merger_Coll_ColorRad3 *= Const_kpc / UNIT_L;

// (2) load the radial profiles
   if ( OPT__INIT != INIT_BY_RESTART ) {

      const std::string filename1(Merger_File_Prof1);
      const std::string filename2(Merger_File_Prof2);
      const std::string filename3(Merger_File_Prof3);

      // cluster 1
      if ( Merger_Coll_IsGas1 ) {

         if ( MPI_Rank == 0 ) {
            Merger_NBin1 = Read_Num_Points_ClusterMerger(filename1);
            Aux_Message(stdout, "num_points1 = %d\n", Merger_NBin1);
         }

         MPI_Bcast(&Merger_NBin1, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R1 = new double [Merger_NBin1];
         Table_D1 = new double [Merger_NBin1];
         Table_P1 = new double [Merger_NBin1];
         Table_M1 = new double [Merger_NBin1];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename1, "/fields/radius", Table_R1);
            Read_Profile_ClusterMerger(filename1, "/fields/density", Table_D1);
            Read_Profile_ClusterMerger(filename1, "/fields/pressure", Table_P1);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename1, "/fields/metallicity", Table_M1);
            else
               for ( int i; i < Merger_NBin1; i++ ) Table_M1[i] = 0.0;

            // convert to code units (assuming the input units are cgs)
            for ( int b=0; b<Merger_NBin1; b++ ) {
               Table_R1[b] /= UNIT_L;
               Table_D1[b] /= UNIT_D;
               Table_P1[b] /= UNIT_P;
            }

         }

         MPI_Bcast(Table_R1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_IsGas1 )

      // cluster 2
      if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2) {

         if (MPI_Rank == 0) {
            Merger_NBin2 = Read_Num_Points_ClusterMerger(filename2);
            Aux_Message(stdout, "num_points2 = %d\n", Merger_NBin2);
         }

         MPI_Bcast(&Merger_NBin2, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R2 = new double [Merger_NBin2];
         Table_D2 = new double [Merger_NBin2];
         Table_P2 = new double [Merger_NBin2];
         Table_M2 = new double [Merger_NBin2];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename2, "/fields/radius", Table_R2);
            Read_Profile_ClusterMerger(filename2, "/fields/density", Table_D2);
            Read_Profile_ClusterMerger(filename2, "/fields/pressure", Table_P2);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename2, "/fields/metallicity", Table_M2);
            else
               for ( int i; i < Merger_NBin2; i++ ) Table_M2[i] = 0.0;
            // convert to code units (assuming the input units are cgs)
            for ( int b=0; b<Merger_NBin2; b++ ) {
               Table_R2[b] /= UNIT_L;
               Table_D2[b] /= UNIT_D;
               Table_P2[b] /= UNIT_P;
            }

         }

         MPI_Bcast(Table_R2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 )

      // cluster 3
      if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3) {

         if ( MPI_Rank == 0 ) {
            Merger_NBin3 = Read_Num_Points_ClusterMerger(filename3);
            Aux_Message(stdout, "num_points3 = %d\n", Merger_NBin3);
         }

         MPI_Bcast(&Merger_NBin3, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R3 = new double [Merger_NBin3];
         Table_D3 = new double [Merger_NBin3];
         Table_P3 = new double [Merger_NBin3];
         Table_M3 = new double [Merger_NBin3];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename3, "/fields/radius", Table_R3);
            Read_Profile_ClusterMerger(filename3, "/fields/density", Table_D3);
            Read_Profile_ClusterMerger(filename3, "/fields/pressure", Table_P3);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename3, "/fields/metallicity", Table_M3);
            else
               for ( int i; i < Merger_NBin3; i++ ) Table_M3[i] = 0.0;
            // convert to code units (assuming the input units are cgs)
            for (int b=0; b<Merger_NBin3; b++) {
               Table_R3[b] /= UNIT_L;
               Table_D3[b] /= UNIT_D;
               Table_P3[b] /= UNIT_P;
            }

         }

         MPI_Bcast(Table_R3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 )

   } // if ( OPT__INIT != INIT_BY_RESTART )


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID        = %d\n",           TESTPROB_ID );
      Aux_Message( stdout, "  number of clusters     = %d\n",           Merger_Coll_NumHalos );
      if ( Merger_Coll_IsGas1 )
      Aux_Message( stdout, "  profile file 1         = %s\n",           Merger_File_Prof1 );
      Aux_Message( stdout, "  particle file 1        = %s\n",           Merger_File_Par1 );
      Aux_Message( stdout, "  cluster 1 w/ gas       = %s\n",          (Merger_Coll_IsGas1)? "yes":"no" );
      if ( Merger_Coll_IsGas1 )
      Aux_Message( stdout, "  cluster 1 color radius = %g\n",           Merger_Coll_ColorRad1 );
      Aux_Message( stdout, "  cluster 1 x-position   = %g\n",           Merger_Coll_PosX1 );
      Aux_Message( stdout, "  cluster 1 y-position   = %g\n",           Merger_Coll_PosY1 );
      Aux_Message( stdout, "  cluster 1 x-velocity   = %g\n",           Merger_Coll_VelX1 );
      Aux_Message( stdout, "  cluster 1 y-velocity   = %g\n",           Merger_Coll_VelY1 );
      if ( Merger_Coll_NumHalos > 1 ) {
      if ( Merger_Coll_IsGas2 )
      Aux_Message( stdout, "  profile file 2         = %s\n",           Merger_File_Prof2 );
      Aux_Message( stdout, "  particle file 2        = %s\n",           Merger_File_Par2 );
      Aux_Message( stdout, "  cluster 2 w/ gas       = %s\n",          (Merger_Coll_IsGas2)? "yes":"no" );
      if ( Merger_Coll_IsGas2 )
      Aux_Message( stdout, "  cluster 2 color radius = %g\n",           Merger_Coll_ColorRad2 );
      Aux_Message( stdout, "  cluster 2 x-position   = %g\n",           Merger_Coll_PosX2 );
      Aux_Message( stdout, "  cluster 2 y-position   = %g\n",           Merger_Coll_PosY2 );
      Aux_Message( stdout, "  cluster 2 x-velocity   = %g\n",           Merger_Coll_VelX2 );
      Aux_Message( stdout, "  cluster 2 y-velocity   = %g\n",           Merger_Coll_VelY2 );
      }
      if ( Merger_Coll_NumHalos > 2 ) {
      if ( Merger_Coll_IsGas3 )
      Aux_Message( stdout, "  profile file 3         = %s\n",           Merger_File_Prof3 );
      Aux_Message( stdout, "  particle file 3        = %s\n",           Merger_File_Par3 );
      Aux_Message( stdout, "  cluster 3 w/ gas       = %s\n",          (Merger_Coll_IsGas3)? "yes":"no" );
      if ( Merger_Coll_IsGas3 )
      Aux_Message( stdout, "  cluster 3 color radius = %g\n",           Merger_Coll_ColorRad3 );
      Aux_Message( stdout, "  cluster 3 x-position   = %g\n",           Merger_Coll_PosX3 );
      Aux_Message( stdout, "  cluster 3 y-position   = %g\n",           Merger_Coll_PosY3 );
      Aux_Message( stdout, "  cluster 3 x-velocity   = %g\n",           Merger_Coll_VelX3 );
      Aux_Message( stdout, "  cluster 3 y-velocity   = %g\n",           Merger_Coll_VelY3 );
      }
      Aux_Message( stdout, "  use metals             = %s\n",          (Merger_Coll_UseMetals)? "yes":"no" );
      Aux_Message( stdout, "=============================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

#endif // #ifdef SUPPORT_HDF5

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

   if ( !( Merger_Coll_IsGas1 || Merger_Coll_IsGas2 || Merger_Coll_IsGas3 ))
      return;

   const double ClusterCenter1[3] = { Merger_Coll_PosX1, Merger_Coll_PosY1, amr->BoxCenter[2] };
   const double ClusterCenter2[3] = { Merger_Coll_PosX2, Merger_Coll_PosY2, amr->BoxCenter[2] };
   const double ClusterCenter3[3] = { Merger_Coll_PosX3, Merger_Coll_PosY3, amr->BoxCenter[2] };

   double r1, r2, r3, Dens1, Dens2, Dens3, Pres1, Pres2, Pres3;
   double Metl1, Metl2, Metl3;
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot, Metl;

   //    for each cell, we sum up the density and pressure from each halo and then calculate the weighted velocity
   if ( Merger_Coll_IsGas1 ) {
      r1    = sqrt( SQR(x-ClusterCenter1[0]) + SQR(y-ClusterCenter1[1]) + SQR(z-ClusterCenter1[2]) );
      Dens1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_D1, r1 );
      Pres1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_P1, r1 );
      Metl1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_M1, r1 );
   } else {
      r1    = HUGE_NUMBER;
      Dens1 = 0.0;
      Pres1 = 0.0;
      Metl1 = 0.0;
   }

   if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 ) {
      r2    = sqrt( SQR(x-ClusterCenter2[0]) + SQR(y-ClusterCenter2[1]) + SQR(z-ClusterCenter2[2]) );
      Dens2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_D2, r2 );
      Pres2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_P2, r2 );
      Metl2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_M2, r2 );
   } else {
      r2    = HUGE_NUMBER;
      Dens2 = 0.0;
      Pres2 = 0.0;
      Metl2 = 0.0;
   }

   if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 ) {
      r3    = sqrt( SQR(x-ClusterCenter3[0]) + SQR(y-ClusterCenter3[1]) + SQR(z-ClusterCenter3[2]) );
      Dens3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_D3, r3 );
      Pres3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_P3, r3 );
      Metl3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_M3, r3 );
   } else {
      r3    = HUGE_NUMBER;
      Dens3 = 0.0;
      Pres3 = 0.0;
      Metl3 = 0.0;
   }

   if ( Dens1 == NULL_REAL )
      Dens1 = Table_D1[Merger_NBin1-1];
   if ( Pres1 == NULL_REAL )
      Pres1 = Table_P1[Merger_NBin1-1];
   if ( Metl1 == NULL_REAL )
      Metl1 = Table_M1[Merger_NBin1-1];

   if ( Dens2 == NULL_REAL )
      Dens2 = Table_D2[Merger_NBin2-1];
   if ( Pres2 == NULL_REAL )
      Pres2 = Table_P2[Merger_NBin2-1];
   if ( Metl2 == NULL_REAL )
      Metl2 = Table_M2[Merger_NBin2-1];

   if ( Dens3 == NULL_REAL )
      Dens3 = Table_D3[Merger_NBin3-1];
   if ( Pres3 == NULL_REAL )
      Pres3 = Table_P3[Merger_NBin3-1];
   if ( Metl3 == NULL_REAL )
      Metl3 = Table_M3[Merger_NBin3-1];

   Dens = Dens1 + Dens2 + Dens3;
   Pres = Pres1 + Pres2 + Pres3;

   MomX = Merger_Coll_VelX1*Dens1 + Merger_Coll_VelX2*Dens2 + Merger_Coll_VelX3*Dens3;
   MomY = Merger_Coll_VelY1*Dens1 + Merger_Coll_VelY2*Dens2 + Merger_Coll_VelY3*Dens3;
   MomZ = 0.0;

   // compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table, NULL );    // assuming EoS requires no passive scalars

   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

   if ( Merger_Coll_UseMetals ) {
      Metl = Metl1*Dens1 + Metl2*Dens2 + Metl3*Dens3;
      fluid[Idx_Metal] = Metl;
   }

   if ( Merger_Coll_IsGas1 ) {
      if ( r1 < Merger_Coll_ColorRad1 )
         fluid[ColorField1Idx] = Dens;
      else
         fluid[ColorField1Idx] = 0.0;
   }

   if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 ) {
      if ( r2 < Merger_Coll_ColorRad2 )
         fluid[ColorField2Idx] = Dens;
      else
         fluid[ColorField2Idx] = 0.0;
   }

   if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 ) {
      if ( r3 < Merger_Coll_ColorRad3 )
         fluid[ColorField3Idx] = Dens;
      else
         fluid[ColorField3Idx] = 0.0;
   }

} // FUNCTION : SetGridIC

#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


//-------------------------------------------------------------------------------------------------------
// Function    :  End_ClusterMerger
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_ClusterMerger()
{

   delete [] Table_R1;
   delete [] Table_D1;
   delete [] Table_P1;
   delete [] Table_M1;

   delete [] Table_R2;
   delete [] Table_D2;
   delete [] Table_P2;
   delete [] Table_M2;

   delete [] Table_R3;
   delete [] Table_D3;
   delete [] Table_P3;
   delete [] Table_M3;

} // FUNCTION : End_ClusterMerger

#ifdef MHD

void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   // Setting the B-field in this problem is always handled by reading from a uniform
   // grid in the file B_IC

   return;

} // FUNCTION : SetBFieldIC

#endif // #ifdef MHD


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ClusterMerger
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ClusterMerger()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr         = SetGridIC;
   End_User_Ptr                   = End_ClusterMerger;
   Par_Init_ByFunction_Ptr        = Par_Init_ByFunction_ClusterMerger;
   Init_Field_User_Ptr            = AddNewField_ClusterMerger;
   Par_Init_Attribute_User_Ptr    = AddNewParticleAttribute_ClusterMerger;
#  ifdef MHD
   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif
#  endif // if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ClusterMerger

#ifdef SUPPORT_HDF5

int Read_Num_Points_ClusterMerger(std::string filename)
{

   hid_t   file_id, dataset, dataspace;
   herr_t  status;
   hsize_t dims[1], maxdims[1];

   int rank;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset   = H5Dopen(file_id, "/fields/radius", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   H5Fclose(file_id);

   return (int)dims[0];

} // FUNCTION : Read_Num_Points_ClusterMerger

void Read_Profile_ClusterMerger(std::string filename, std::string fieldname,
                                double field[])
{

   hid_t   file_id, dataset;
   herr_t  status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset = H5Dopen(file_id, fieldname.c_str(), H5P_DEFAULT);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, field);
   H5Dclose(dataset);

   H5Fclose(file_id);

   return;

} // FUNCTION : Read_Profile_ClusterMerger

#endif // #ifdef SUPPORT_HDF5

void AddNewField_ClusterMerger()
{

   if ( Merger_Coll_UseMetals )
      Idx_Metal = AddField( "Metal", NORMALIZE_NO );
   if ( ColorField1Idx == Idx_Undefined )
      ColorField1Idx = AddField( "ColorField1", NORMALIZE_NO );
   if ( Merger_Coll_NumHalos > 1 && ColorField2Idx == Idx_Undefined )
      ColorField2Idx = AddField( "ColorField2", NORMALIZE_NO );
   if ( Merger_Coll_NumHalos > 2 && ColorField3Idx == Idx_Undefined )
      ColorField3Idx = AddField( "ColorField3", NORMALIZE_NO );

}
