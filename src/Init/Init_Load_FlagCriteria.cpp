#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_FlagCriteria
// Description :  Load the flag criteria from several input files named "Input__Flag_XXX"
//
// Note        :  The loaded data will be used as the refinement thresholds at each refinement level
//-------------------------------------------------------------------------------------------------------
void Init_Load_FlagCriteria()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  if ( MODEL != HYDRO )
   const bool OPT__FLAG_PRES_GRADIENT = false;
   const bool OPT__FLAG_ENGY_GRADIENT = false;
   double *FlagTable_PresGradient     = NULL;
   double *FlagTable_EngyGradient     = NULL;

   const bool OPT__FLAG_VORTICITY     = false;
   double *FlagTable_Vorticity        = NULL;

   const bool OPT__FLAG_JEANS         = false;
   double *FlagTable_Jeans            = NULL;
#  endif

#  ifndef MHD
   const bool OPT__FLAG_CURRENT       = false;
   double *FlagTable_Current          = NULL;
#  endif

#  ifndef SRHD
   const bool OPT__FLAG_LORENTZ_GRADIENT   = false;
   double *FlagTable_LorentzFactorGradient = NULL;


   const bool OPT__FLAG_4VELOCITY          = false;
   double *FlagTable_4Velocity             = NULL;

   const bool OPT__FLAG_MOM_OVER_DENS      = false;
   double *FlagTable_Mom_Over_Dens         = NULL;
#  endif


#  if ( MODEL != ELBDM )
   const bool OPT__FLAG_ENGY_DENSITY  = false;
   double FlagTable_EngyDensity[NLEVEL-1][2];
#  endif

#  ifndef PARTICLE
   const bool OPT__FLAG_NPAR_PATCH    = false;
   int *FlagTable_NParPatch           = NULL;

   const bool OPT__FLAG_NPAR_CELL     = false;
   int *FlagTable_NParCell            = NULL;

   const bool OPT__FLAG_PAR_MASS_CELL = false;
   double *FlagTable_ParMassCell      = NULL;
#  endif

#  if   ( MODEL == HYDRO )
   const bool OPT__FLAG_LOHNER = ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP );
#  elif ( MODEL == ELBDM )
   const bool OPT__FLAG_LOHNER = OPT__FLAG_LOHNER_DENS;
#  else
#  error : unsupported MODEL !!
#  endif

   const int  NFlagMode         = 16;
   const bool Flag[NFlagMode]   = { OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_PRES_GRADIENT,
                                    OPT__FLAG_ENGY_DENSITY, OPT__FLAG_LOHNER, OPT__FLAG_USER,
                                    (bool)OPT__FLAG_NPAR_PATCH, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL,
                                    OPT__FLAG_VORTICITY, OPT__FLAG_JEANS, OPT__FLAG_CURRENT, OPT__FLAG_ENGY_GRADIENT,
                                    OPT__FLAG_4VELOCITY, OPT__FLAG_MOM_OVER_DENS, OPT__FLAG_LORENTZ_GRADIENT };
   const char ModeName[][100]   = { "OPT__FLAG_RHO", "OPT__FLAG_RHO_GRADIENT", "OPT__FLAG_PRES_GRADIENT",
                                    "OPT__FLAG_ENGY_DENSITY", "OPT__FLAG_LOHNER", "OPT__FLAG_USER",
                                    "OPT__FLAG_NPAR_PATCH", "OPT__FLAG_NPAR_CELL", "OPT__FLAG_PAR_MASS_CELL",
                                    "OPT__FLAG_VORTICITY", "OPT__FLAG_JEANS", "OPT__FLAG_CURRENT", "OPT__FLAG_ENGY_GRADIENT",
                                    "OPT__FLAG_4VELOCITY", "OPT__FLAG_MOM_OVER_DENS", "OPT__FLAG_LORENTZ_GRADIENT" };
   const char FileName[][100]   = { "Input__Flag_Rho", "Input__Flag_RhoGradient", "Input__Flag_PresGradient",
                                    "Input__Flag_EngyDensity", "Input__Flag_Lohner", "Input__Flag_User",
                                    "Input__Flag_NParPatch", "Input__Flag_NParCell", "Input__Flag_ParMassCell",
                                    "Input__Flag_Vorticity", "Input__Flag_Jeans", "Input__Flag_Current", "Input__Flag_EngyGradient",
                                    "Input__Flag_4Velocity", "Input__Flag_Mom_Over_Dens", "Input__Flag_LorentzGradient" };
   double *FlagTable[NFlagMode] = { FlagTable_Rho, FlagTable_RhoGradient, FlagTable_PresGradient,
                                    NULL, NULL, NULL, NULL, NULL, FlagTable_ParMassCell,
                                    FlagTable_Vorticity, FlagTable_Jeans, FlagTable_Current, FlagTable_EngyGradient,
                                    FlagTable_4Velocity, FlagTable_Mom_Over_Dens, FlagTable_LorentzFactorGradient };

   FILE *File;
   char *input_line = NULL, TargetName[100];
   size_t len = 0;
   int Trash, n;


// initialize the flag tables as -1
   for (int lv=0; lv<NLEVEL-1; lv++)
   {
      FlagTable_Rho         [lv]    = -1.0;
      FlagTable_RhoGradient [lv]    = -1.0;

      for (int t=0; t<4; t++)
      FlagTable_Lohner      [lv][t] = -1.0;

      FlagTable_User        [lv]    = (double *) malloc( OPT__FLAG_USER_NUM*sizeof(double) );
      for (int t=0; t<OPT__FLAG_USER_NUM; t++)
      FlagTable_User        [lv][t] = -1.0;

#     if   ( MODEL == HYDRO )
      FlagTable_PresGradient[lv]    = -1.0;
      FlagTable_EngyGradient[lv]    = -1.0;
      FlagTable_Vorticity   [lv]    = -1.0;
      FlagTable_Jeans       [lv]    = -1.0;
#     ifdef MHD
      FlagTable_Current     [lv]    = -1.0;
#     endif
#     ifdef SRHD
      FlagTable_4Velocity             [lv]    = -1.0;
      FlagTable_Mom_Over_Dens         [lv]    = -1.0;
      FlagTable_LorentzFactorGradient [lv]    = -1.0;
#     endif

#     elif ( MODEL == ELBDM )
      for (int t=0; t<2; t++)
      FlagTable_EngyDensity [lv][t] = -1.0;
#     endif

#     ifdef PARTICLE
      FlagTable_NParPatch   [lv]    = -1;
      FlagTable_NParCell    [lv]    = -1;
      FlagTable_ParMassCell [lv]    = -1.0;
#     endif
   } // for (int lv=0; lv<NLEVEL-1; lv++)


// nothing to do if there is no refinement level
   if ( MAX_LEVEL > 0 )
   for (int FlagMode=0; FlagMode<NFlagMode; FlagMode++)
   {
      if ( Flag[FlagMode] )
      {
//       open the target file
         strcpy( TargetName, FileName[FlagMode] );

         if ( !Aux_CheckFileExist(TargetName) )
         {
            Aux_Message( stderr, "\n" );
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist for the mode \"%s\" !!\n",
                       TargetName, ModeName[FlagMode] );
         }

         File = fopen( TargetName, "r" );

//       skip the header
         getline( &input_line, &len, File );

//       begin to read
         for (int lv=0; lv<MAX_LEVEL; lv++)
         {
            n = getline( &input_line, &len, File );

//          check
            if ( n <= 1 )
            {
               Aux_Message( stderr, "\n" );
               Aux_Error( ERROR_INFO, "incorrect reading at level %d of the file <%s> !!\n",
                          lv, TargetName );
            }

//          OPT__FLAG_ENGY_DENSITY and OPT__FLAG_LOHNER have two and three columns to be loaded, respectively
            if      ( FlagMode == 3 )  sscanf( input_line, "%d%lf%lf", &Trash, &FlagTable_EngyDensity[lv][0],
                                                                               &FlagTable_EngyDensity[lv][1] );
            else if ( FlagMode == 4 )  sscanf( input_line, "%d%lf%lf%lf%lf", &Trash, &FlagTable_Lohner[lv][0],
                                                                                     &FlagTable_Lohner[lv][1],
                                                                                     &FlagTable_Lohner[lv][2],
                                                                                     &FlagTable_Lohner[lv][3] );
//          OPT__FLAG_USER has OPT__FLAG_USER_NUM columns to be loaded
            else if ( FlagMode == 5 )
            {
               char *input_line0 = input_line;
               int offset;

               sscanf( input_line, "%d%n", &Trash, &offset );
               input_line += offset;

               for (int t=0; t<OPT__FLAG_USER_NUM; t++)
               {
                  sscanf( input_line, "%lf%n", &FlagTable_User[lv][t], &offset );
                  input_line += offset;
               }

               input_line = input_line0;  // reset input_line
            }
//          OPT__FLAG_NPAR_PATCH/CELL load integers
            else if ( FlagMode == 6 )  sscanf( input_line, "%d%d",  &Trash, &FlagTable_NParPatch[lv] );
            else if ( FlagMode == 7 )  sscanf( input_line, "%d%d",  &Trash, &FlagTable_NParCell [lv] );

//          others use the default format: (integer, double)
            else                       sscanf( input_line, "%d%lf", &Trash, &FlagTable[FlagMode][lv] );
         }

         fclose( File );

      } // if ( Flag[FlagMode] )
   } // for (int FlagMode=0; FlagMode<NFlagMode; FlagMode++)


   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Load_FlagCriteria
