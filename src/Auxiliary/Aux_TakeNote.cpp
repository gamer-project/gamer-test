#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif
#include <sched.h>
#ifdef __APPLE__
#include <cpuid.h>
#endif

static int get_cpuid();




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_TakeNote
// Description :  Record simulation parameters and the content in the file "Input__Note" to the
//                note file "Record__Note"
//-------------------------------------------------------------------------------------------------------
void Aux_TakeNote()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_TakeNote ...\n" );


   const char FileName[] = "Record__Note";
   FILE *Note;

   if ( MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

//    copy the content in the file "Input__Note"
      Note = fopen( FileName, "a" );
      fprintf( Note, "\n\n\nSimulation Notes\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fclose( Note );

      system( "cat ./Input__Note >> Record__Note" );

      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in the Makefile (numerical schemes)
      fprintf( Note, "Makefile Options (numerical schemes)\n" );
      fprintf( Note, "***********************************************************************************\n" );

//    a. options for all models
#     if   ( MODEL == HYDRO )
      fprintf( Note, "MODEL                           HYDRO\n" );
#     elif ( MODEL == ELBDM )
      fprintf( Note, "MODEL                           ELBDM\n" );
#     elif ( MODEL == PAR_ONLY )
      fprintf( Note, "MODEL                           PAR_ONLY\n" );
#     else
      fprintf( Note, "MODEL                           UNKNOWN\n" );
#     endif // MODEL

#     ifdef GRAVITY
      fprintf( Note, "GRAVITY                         ON\n" );
#     else
      fprintf( Note, "GRAVITY                         OFF\n" );
#     endif

#     ifdef GRAVITY
#     if   ( POT_SCHEME == SOR )
      fprintf( Note, "POT_SCHEME                      SOR\n" );
#     elif ( POT_SCHEME == MG )
      fprintf( Note, "POT_SCHEME                      MG\n" );
#     elif ( POT_SCHEME == NONE )
      fprintf( Note, "POT_SCHEME                      NONE\n" );
#     else
      fprintf( Note, "POT_SCHEME                      UNKNOWN\n" );
#     endif

#     ifdef STORE_POT_GHOST
      fprintf( Note, "STORE_POT_GHOST                 ON\n" );
#     else
      fprintf( Note, "STORE_POT_GHOST                 OFF\n" );
#     endif

#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "UNSPLIT_GRAVITY                 ON\n" );
#     else
      fprintf( Note, "UNSPLIT_GRAVITY                 OFF\n" );
#     endif
#     endif // #ifdef GRAVITY

#     ifdef COMOVING
      fprintf( Note, "COMOVING                        ON\n" );
#     else
      fprintf( Note, "COMOVING                        OFF\n" );
#     endif

#     ifdef PARTICLE
      fprintf( Note, "PARTICLE                        ON\n" );
#     else
      fprintf( Note, "PARTICLE                        OFF\n" );
#     endif

#     ifdef SUPPORT_GRACKLE
      fprintf( Note, "SUPPORT_GRACKLE                 ON\n" );
#     else
      fprintf( Note, "SUPPORT_GRACKLE                 OFF\n" );
#     endif

//    b. options in HYDRO
#     if   ( MODEL == HYDRO )

#     if   ( FLU_SCHEME == RTVD )
      fprintf( Note, "FLU_SCHEME                      RTVD\n" );
#     elif ( FLU_SCHEME == MHM )
      fprintf( Note, "FLU_SCHEME                      MHM\n" );
#     elif ( FLU_SCHEME == MHM_RP )
      fprintf( Note, "FLU_SCHEME                      MHM with Riemann prediction\n" );
#     elif ( FLU_SCHEME == CTU )
      fprintf( Note, "FLU_SCHEME                      CTU\n" );
#     elif ( FLU_SCHEME == NONE )
      fprintf( Note, "FLU_SCHEME                      NONE\n" );
#     else
      fprintf( Note, "FLU_SCHEME                      UNKNOWN\n" );
#     endif

#     if   ( LR_SCHEME == PLM )
      fprintf( Note, "LR_SCHEME                       PLM\n" );
#     elif ( LR_SCHEME == PPM )
      fprintf( Note, "LR_SCHEME                       PPM\n" );
#     elif ( LR_SCHEME == NONE )
      fprintf( Note, "LR_SCHEME                       NONE\n" );
#     else
      fprintf( Note, "LR_SCHEME                       UNKNOWN\n" );
#     endif

#     if   ( RSOLVER == EXACT )
      fprintf( Note, "RSOLVER                         EXACT\n" );
#     elif ( RSOLVER == ROE )
      fprintf( Note, "RSOLVER                         ROE\n" );
#     elif ( RSOLVER == HLLE )
      fprintf( Note, "RSOLVER                         HLLE\n" );
#     elif ( RSOLVER == HLLC )
      fprintf( Note, "RSOLVER                         HLLC\n" );
#     elif ( RSOLVER == HLLD )
      fprintf( Note, "RSOLVER                         HLLD\n" );
#     elif ( RSOLVER == NONE )
      fprintf( Note, "RSOLVER                         NONE\n" );
#     else
      fprintf( Note, "RSOLVER                         UNKNOWN\n" );
#     endif

#     if   ( DUAL_ENERGY == DE_ENPY )
      fprintf( Note, "DUAL_ENERGY                     DE_ENPY\n" );
#     elif ( DUAL_ENERGY == DE_EINT )
      fprintf( Note, "DUAL_ENERGY                     DE_EINT\n" );
#     elif ( DUAL_ENERGY == NONE )
      fprintf( Note, "DUAL_ENERGY                     NONE\n" );
#     else
      fprintf( Note, "DUAL_ENERGY                     UNKNOWN\n" );
#     endif

#     ifdef MHD
      fprintf( Note, "MHD                             ON\n" );
#     else
      fprintf( Note, "MHD                             OFF\n" );
#     endif

#     ifdef SRHD
      fprintf( Note, "SRHD                            ON\n" );
#     ifdef REDUCED_ENERGY
      fprintf( Note, "REDUCED_ENERGY                  ON\n" );
#     else
      fprintf( Note, "REDUCED_ENERGY                  OFF\n" );
#     endif

#     ifdef FOUR_VELOCITY
      fprintf( Note, "FOUR_VELOCITY                  ON\n" );
#     else
      fprintf( Note, "FOUR_VELOCITY                  OFF\n" );
#     endif
#     else
      fprintf( Note, "SRHD                            OFF\n" );
#     endif

#     ifdef COSMIC_RAY
      fprintf( Note, "COSMIC_RAY                      ON\n" );
#     else
      fprintf( Note, "COSMIC_RAY                      OFF\n" );
#     endif

#     if   ( EOS == EOS_GAMMA )
      fprintf( Note, "EOS                             EOS_GAMMA\n" );
#     elif ( EOS == EOS_ISOTHERMAL )
      fprintf( Note, "EOS                             EOS_ISOTHERMAL\n" );
#     elif ( EOS == EOS_NUCLEAR )
      fprintf( Note, "EOS                             EOS_NUCLEAR\n" );
#     elif ( EOS == EOS_TAUBMATHEWS )
      fprintf( Note, "EOS                             EOS_TAUBMATHEWS\n" );
#     elif ( EOS == EOS_TABULAR )
      fprintf( Note, "EOS                             EOS_TABULAR\n" );
#     elif ( EOS == EOS_USER )
      fprintf( Note, "EOS                             EOS_USER\n" );
#     else
      fprintf( Note, "EOS                             UNKNOWN\n" );
#     endif

#     ifdef BAROTROPIC_EOS
      fprintf( Note, "BAROTROPIC_EOS                  ON\n" );
#     else
      fprintf( Note, "BAROTROPIC_EOS                  OFF\n" );
#     endif

//    c. options in ELBDM
#     elif ( MODEL == ELBDM )

#     ifdef CONSERVE_MASS
      fprintf( Note, "CONSERVE_MASS                   ON\n" );
#     else
      fprintf( Note, "CONSERVE_MASS                   OFF\n" );
#     endif

#     ifdef LAPLACIAN_4TH
      fprintf( Note, "LAPLACIAN_4TH                   ON\n" );
#     else
      fprintf( Note, "LAPLACIAN_4TH                   OFF\n" );
#     endif

#     ifdef QUARTIC_SELF_INTERACTION
      fprintf( Note, "QUARTIC_SELF_INTERACTION        ON\n" );
#     else
      fprintf( Note, "QUARTIC_SELF_INTERACTION        OFF\n" );
#     endif

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

//    d. options in PARTICLE
#     ifdef PARTICLE
#     ifdef STORE_PAR_ACC
      fprintf( Note, "STORE_PAR_ACC                   ON\n" );
#     else
      fprintf( Note, "STORE_PAR_ACC                   OFF\n" );
#     endif

#     ifdef STAR_FORMATION
      fprintf( Note, "STAR_FORMATION                  ON\n" );
#     else
      fprintf( Note, "STAR_FORMATION                  OFF\n" );
#     endif
#     endif

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in the Makefile (optimization and compilation)
      fprintf( Note, "Makefile Options (optimization and compilation)\n" );
      fprintf( Note, "***********************************************************************************\n" );

#     ifdef GPU
      fprintf( Note, "GPU                             ON\n" );
#     else
      fprintf( Note, "GPU                             OFF\n" );
#     endif

#     ifdef GAMER_DEBUG
      fprintf( Note, "GAMER_DEBUG                     ON\n" );
#     else
      fprintf( Note, "GAMER_DEBUG                     OFF\n" );
#     endif

#     ifdef BITWISE_REPRODUCIBILITY
      fprintf( Note, "BITWISE_REPRODUCIBILITY         ON\n" );
#     else
      fprintf( Note, "BITWISE_REPRODUCIBILITY         OFF\n" );
#     endif

#     ifdef TIMING
      fprintf( Note, "TIMING                          ON\n" );
#     else
      fprintf( Note, "TIMING                          OFF\n" );
#     endif

#     ifdef TIMING_SOLVER
      fprintf( Note, "TIMING_SOLVER                   ON\n" );
#     else
      fprintf( Note, "TIMING_SOLVER                   OFF\n" );
#     endif

#     ifdef FLOAT8
      fprintf( Note, "FLOAT8                          ON\n" );
#     else
      fprintf( Note, "FLOAT8                          OFF\n" );
#     endif

#     ifdef SERIAL
      fprintf( Note, "SERIAL                          ON\n" );
#     else
      fprintf( Note, "SERIAL                          OFF\n" );
#     endif

#     ifdef LOAD_BALANCE
#     if   ( LOAD_BALANCE == HILBERT )
      fprintf( Note, "LOAD_BALANCE                    HILBERT\n" );
#     else
      fprintf( Note, "LOAD_BALANCE                    UNKNOWN\n" );
#     endif
#     else // #ifdef LOAD_BALANCE
      fprintf( Note, "LOAD_BALANCE                    OFF\n" );
#     endif // #ifdef LOAD_BALANCE ... else ...

#     ifdef OVERLAP_MPI
      fprintf( Note, "OVERLAP_MPI                     ON\n" );
#     else
      fprintf( Note, "OVERLAP_MPI                     OFF\n" );
#     endif

#     ifdef OPENMP
      fprintf( Note, "OPENMP                          ON\n" );
#     else
      fprintf( Note, "OPENMP                          OFF\n" );
#     endif

#     ifdef GPU
#     if   ( GPU_ARCH == FERMI )
      fprintf( Note, "GPU_ARCH                        FERMI\n" );
#     elif ( GPU_ARCH == KEPLER )
      fprintf( Note, "GPU_ARCH                        KEPLER\n" );
#     elif ( GPU_ARCH == MAXWELL )
      fprintf( Note, "GPU_ARCH                        MAXWELL\n" );
#     elif ( GPU_ARCH == PASCAL )
      fprintf( Note, "GPU_ARCH                        PASCAL\n" );
#     elif ( GPU_ARCH == VOLTA )
      fprintf( Note, "GPU_ARCH                        VOLTA\n" );
#     elif ( GPU_ARCH == TURING )
      fprintf( Note, "GPU_ARCH                        TURING\n" );
#     else
      fprintf( Note, "GPU_ARCH                        UNKNOWN\n" );
#     endif
#     endif

#     ifdef LAOHU
      fprintf( Note, "LAOHU                           ON\n" );
#     else
      fprintf( Note, "LAOHU                           OFF\n" );
#     endif

#     ifdef SUPPORT_HDF5
      fprintf( Note, "SUPPORT_HDF5                    ON\n" );
#     else
      fprintf( Note, "SUPPORT_HDF5                    OFF\n" );
#     endif

#     ifdef SUPPORT_GSL
      fprintf( Note, "SUPPORT_GSL                     ON\n" );
#     else
      fprintf( Note, "SUPPORT_GSL                     OFF\n" );
#     endif

#     if   ( RANDOM_NUMBER == RNG_GNU_EXT )
      fprintf( Note, "RANDOM_NUMBER                   RNG_GNU_EXT\n" );
#     elif ( RANDOM_NUMBER == RNG_CPP11 )
      fprintf( Note, "RANDOM_NUMBER                   RNG_CPP11\n" );
#     else
      fprintf( Note, "RANDOM_NUMBER                   UNKNOWN\n" );
#     endif

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the simulation options in CUFLU.h and CUPOT.h
      fprintf( Note, "Other Options (in CUFLU.h and CUPOT.h)\n" );
      fprintf( Note, "***********************************************************************************\n" );

#     ifdef BIT_REP_FLUX
      fprintf( Note, "BIT_REP_FLUX                    ON\n" );
#     else
      fprintf( Note, "BIT_REP_FLUX                    OFF\n" );
#     endif

#     ifdef MHD
#     ifdef BIT_REP_ELECTRIC
      fprintf( Note, "BIT_REP_ELECTRIC                ON\n" );
#     else
      fprintf( Note, "BIT_REP_ELECTRIC                OFF\n" );
#     endif
#     endif

#     if   ( MODEL == HYDRO )
#     ifdef CHECK_NEGATIVE_IN_FLUID
      fprintf( Note, "CHECK_NEGATIVE_IN_FLUID         ON\n" );
#     else
      fprintf( Note, "CHECK_NEGATIVE_IN_FLUID         OFF\n" );
#     endif

#     ifdef CHAR_RECONSTRUCTION
      fprintf( Note, "CHAR_RECONSTRUCTION             ON\n" );
#     else
      fprintf( Note, "CHAR_RECONSTRUCTION             OFF\n" );
#     endif

#     ifdef LR_EINT
      fprintf( Note, "LR_EINT                         ON\n" );
#     else
      fprintf( Note, "LR_EINT                         OFF\n" );
#     endif

#     if   ( CHECK_INTERMEDIATE == EXACT )
      fprintf( Note, "CHECK_INTERMEDIATE              EXACT\n" );
#     elif ( CHECK_INTERMEDIATE == HLLE )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLE\n" );
#     elif ( CHECK_INTERMEDIATE == HLLC )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLC\n" );
#     elif ( CHECK_INTERMEDIATE == HLLD )
      fprintf( Note, "CHECK_INTERMEDIATE              HLLD\n" );
#     elif ( CHECK_INTERMEDIATE == NONE )
      fprintf( Note, "CHECK_INTERMEDIATE              OFF\n" );
#     else
      fprintf( Note, "CHECK_INTERMEDIATE              UNKNOWN\n" );
#     endif

#     ifdef HLL_NO_REF_STATE
      fprintf( Note, "HLL_NO_REF_STATE                ON\n" );
#     else
      fprintf( Note, "HLL_NO_REF_STATE                OFF\n" );
#     endif

#     ifdef HLL_INCLUDE_ALL_WAVES
      fprintf( Note, "HLL_INCLUDE_ALL_WAVES           ON\n" );
#     else
      fprintf( Note, "HLL_INCLUDE_ALL_WAVES           OFF\n" );
#     endif

      fprintf( Note, "HLLC_WAVESPEED                  %d\n",      HLLC_WAVESPEED );
      fprintf( Note, "HLLE_WAVESPEED                  %d\n",      HLLE_WAVESPEED );
#     ifdef MHD
      fprintf( Note, "HLLD_WAVESPEED                  %d\n",      HLLD_WAVESPEED );
#     endif

#     ifdef MHD
#     ifdef EULERY
      fprintf( Note, "EULERY                          ON\n" );
#     else
      fprintf( Note, "EULERY                          OFF\n" );
#     endif
#     endif // #ifdef MHD

#     elif ( MODEL == ELBDM )

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

#     if ( defined GRAVITY  &&  POT_SCHEME == SOR  &&  defined GPU )
#     ifdef USE_PSOLVER_10TO14
      fprintf( Note, "USE_PSOLVER_10TO14              ON\n" );
#     else
      fprintf( Note, "USE_PSOLVER_10TO14              OFF\n" );
#     endif

#     ifdef USE_PSOLVER_10TO14
#     ifdef SOR_RHO_SHARED
      fprintf( Note, "SOR_RHO_SHARED                  ON\n" );
#     else
      fprintf( Note, "SOR_RHO_SHARED                  OFF\n" );
#     endif

#     ifdef SOR_CPOT_SHARED
      fprintf( Note, "SOR_CPOT_SHARED                 ON\n" );
#     else
      fprintf( Note, "SOR_CPOT_SHARED                 OFF\n" );
#     endif

#     ifdef SOR_USE_SHUFFLE
      fprintf( Note, "SOR_USE_SHUFFLE                 ON\n" );
#     else
      fprintf( Note, "SOR_USE_SHUFFLE                 OFF\n" );
#     endif

#     ifdef SOR_USE_PADDING
      fprintf( Note, "SOR_USE_PADDING                 ON\n" );
#     else
      fprintf( Note, "SOR_USE_PADDING                 OFF\n" );
#     endif

      fprintf( Note, "SOR_MOD_REDUCTION               %d\n",      SOR_MOD_REDUCTION       );
#     endif // #ifdef USE_PSOLVER_10TO14
#     endif // #if ( defined GRAVITY  &&  POT_SCHEME == SOR  &&  defined GPU )

#     ifdef GPU
#     ifdef DT_FLU_USE_SHUFFLE
      fprintf( Note, "DT_FLU_USE_SHUFFLE              ON\n" );
#     else
      fprintf( Note, "DT_FLU_USE_SHUFFLE              OFF\n" );
#     endif
#     ifdef GRAVITY
#     ifdef DT_GRA_USE_SHUFFLE
      fprintf( Note, "DT_GRA_USE_SHUFFLE              ON\n" );
#     else
      fprintf( Note, "DT_GRA_USE_SHUFFLE              OFF\n" );
#     endif
#     endif // #ifdef GRAVITY
#     endif // #ifdef GPU

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the symbolic constants
      fprintf( Note, "Symbolic Constants\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "#define VERSION                 %s\n",      VERSION               );
      fprintf( Note, "#define NCOMP_FLUID             %d\n",      NCOMP_FLUID           );
      fprintf( Note, "#define NCOMP_PASSIVE           %d\n",      NCOMP_PASSIVE         );
      fprintf( Note, "#define FLU_NIN                 %d\n",      FLU_NIN               );
      fprintf( Note, "#define FLU_NOUT                %d\n",      FLU_NOUT              );
      fprintf( Note, "#define FLU_NIN_T               %d\n",      FLU_NIN_T             );
      fprintf( Note, "#define FLU_NIN_S               %d\n",      FLU_NIN_S             );
      fprintf( Note, "#define FLU_NOUT_S              %d\n",      FLU_NOUT_S            );
      fprintf( Note, "#define NFLUX_FLUID             %d\n",      NFLUX_FLUID           );
      fprintf( Note, "#define NFLUX_PASSIVE           %d\n",      NFLUX_PASSIVE         );
#     ifdef GRAVITY
      fprintf( Note, "#define GRA_NIN                 %d\n",      GRA_NIN               );
#     endif
#     ifdef MHD
      fprintf( Note, "#define NCOMP_MAG               %d\n",      NCOMP_MAG             );
      fprintf( Note, "#define NCOMP_ELE               %d\n",      NCOMP_ELE             );
#     endif
      fprintf( Note, "#define PATCH_SIZE              %d\n",      PATCH_SIZE            );
      fprintf( Note, "#define MAX_PATCH               %d\n",      MAX_PATCH             );
      fprintf( Note, "#define NLEVEL                  %d\n",      NLEVEL                );
      fprintf( Note, "\n" );
      fprintf( Note, "#define FLU_GHOST_SIZE          %d\n",      FLU_GHOST_SIZE        );
#     if ( MODEL == HYDRO  &&  defined LR_GHOST_SIZE )
      fprintf( Note, "#define LR_GHOST_SIZE           %d\n",      LR_GHOST_SIZE         );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "#define POT_GHOST_SIZE          %d\n",      POT_GHOST_SIZE        );
      fprintf( Note, "#define RHO_GHOST_SIZE          %d\n",      RHO_GHOST_SIZE        );
      fprintf( Note, "#define GRA_GHOST_SIZE          %d\n",      GRA_GHOST_SIZE        );
#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "#define USG_GHOST_SIZE_F        %d\n",      USG_GHOST_SIZE_F      );
      fprintf( Note, "#define USG_GHOST_SIZE_G        %d\n",      USG_GHOST_SIZE_G      );
#     endif
#     ifdef PARTICLE
      fprintf( Note, "#define RHOEXT_GHOST_SIZE       %d\n",      RHOEXT_GHOST_SIZE     );
#     endif
#     endif
      fprintf( Note, "#define SRC_GHOST_SIZE          %d\n",      SRC_GHOST_SIZE        );
      fprintf( Note, "#define FLU_NXT                 %d\n",      FLU_NXT               );
#     ifdef GRAVITY
      fprintf( Note, "#define POT_NXT                 %d\n",      POT_NXT               );
      fprintf( Note, "#define RHO_NXT                 %d\n",      RHO_NXT               );
      fprintf( Note, "#define GRA_NXT                 %d\n",      GRA_NXT               );
#     ifdef UNSPLIT_GRAVITY
      fprintf( Note, "#define USG_NXT_F               %d\n",      USG_NXT_F             );
      fprintf( Note, "#define USG_NXT_G               %d\n",      USG_NXT_G             );
#     endif
#     endif
#     ifdef PARTICLE
      fprintf( Note, "#define RHOEXT_NXT              %d\n",      RHOEXT_NXT            );
#     endif
      fprintf( Note, "#define SRC_NXT                 %d\n",      SRC_NXT               );
#     if ( MODEL == HYDRO )
      fprintf( Note, "#define EOS_NAUX_MAX            %d\n",      EOS_NAUX_MAX          );
      fprintf( Note, "#define EOS_NTABLE_MAX          %d\n",      EOS_NTABLE_MAX        );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "#define EXT_POT_NAUX_MAX        %d\n",      EXT_POT_NAUX_MAX      );
      fprintf( Note, "#define EXT_ACC_NAUX_MAX        %d\n",      EXT_ACC_NAUX_MAX      );
#     endif
      fprintf( Note, "#define SRC_NAUX_DLEP           %d\n",      SRC_NAUX_DLEP         );
      fprintf( Note, "#define SRC_DLEP_PROF_NVAR      %d\n",      SRC_DLEP_PROF_NVAR    );
      fprintf( Note, "#define SRC_DLEP_PROF_NBINMAX   %d\n",      SRC_DLEP_PROF_NBINMAX );
      fprintf( Note, "#define SRC_NAUX_USER           %d\n",      SRC_NAUX_USER         );
#     ifdef GPU
      fprintf( Note, "#define FLU_BLOCK_SIZE_X        %d\n",      FLU_BLOCK_SIZE_X      );
      fprintf( Note, "#define FLU_BLOCK_SIZE_Y        %d\n",      FLU_BLOCK_SIZE_Y      );
#     ifdef GRAVITY
#     if   ( POT_SCHEME == SOR )
      fprintf( Note, "#define POT_BLOCK_SIZE_Z        %d\n",      POT_BLOCK_SIZE_Z      );
#     elif ( POT_SCHEME == MG )
      fprintf( Note, "#define POT_BLOCK_SIZE_X        %d\n",      POT_BLOCK_SIZE_X      );
#     endif
      fprintf( Note, "#define EXTPOT_BLOCK_SIZE       %d\n",      EXTPOT_BLOCK_SIZE     );
      fprintf( Note, "#define GRA_BLOCK_SIZE          %d\n",      GRA_BLOCK_SIZE        );
#     endif // #ifdef GRAVITY
      fprintf( Note, "#define DT_FLU_BLOCK_SIZE       %d\n",      DT_FLU_BLOCK_SIZE     );
#     ifdef GRAVITY
      fprintf( Note, "#define DT_GRA_BLOCK_SIZE       %d\n",      DT_GRA_BLOCK_SIZE     );
#     endif
      fprintf( Note, "#define SRC_BLOCK_SIZE          %d\n",      SRC_BLOCK_SIZE        );
#     endif // #ifdef GPU
#     ifdef PARTICLE
      fprintf( Note, "#define PAR_NATT_TOTAL          %d\n",      PAR_NATT_TOTAL        );
      fprintf( Note, "#define PAR_NATT_USER           %d\n",      PAR_NATT_USER         );
      fprintf( Note, "#define PAR_NATT_STORED         %d\n",      PAR_NATT_STORED       );
#     endif
      fprintf( Note, "#define MAX_STRING              %d\n",      MAX_STRING            );
      fprintf( Note, "#define TINY_NUMBER             %20.14e\n", TINY_NUMBER           );
      fprintf( Note, "#define HUGE_NUMBER             %20.14e\n", HUGE_NUMBER           );
      fprintf( Note, "#define MAX_ERROR               %20.14e\n", MAX_ERROR             );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of simulation scale
      fprintf( Note, "Parameters of Simulation Scale\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "BOX_SIZE (input)                %20.14e\n", BOX_SIZE         );
      fprintf( Note, "BOX_SIZE_X                      %20.14e\n", amr->BoxSize[0]  );
      fprintf( Note, "BOX_SIZE_Y                      %20.14e\n", amr->BoxSize[1]  );
      fprintf( Note, "BOX_SIZE_Z                      %20.14e\n", amr->BoxSize[2]  );
      fprintf( Note, "BOX_SCALE_X                     %d\n",      amr->BoxScale[0] );
      fprintf( Note, "BOX_SCALE_Y                     %d\n",      amr->BoxScale[1] );
      fprintf( Note, "BOX_SCALE_Z                     %d\n",      amr->BoxScale[2] );
      fprintf( Note, "NX0_TOT[0]                      %d\n",      NX0_TOT[0]       );
      fprintf( Note, "NX0_TOT[1]                      %d\n",      NX0_TOT[1]       );
      fprintf( Note, "NX0_TOT[2]                      %d\n",      NX0_TOT[2]       );
      fprintf( Note, "MPI_NRank                       %d\n",      MPI_NRank        );
      fprintf( Note, "MPI_NRank_X[0]                  %d\n",      MPI_NRank_X[0]   );
      fprintf( Note, "MPI_NRank_X[1]                  %d\n",      MPI_NRank_X[1]   );
      fprintf( Note, "MPI_NRank_X[2]                  %d\n",      MPI_NRank_X[2]   );
      fprintf( Note, "OMP_NTHREAD                     %d\n",      OMP_NTHREAD      );
      fprintf( Note, "END_T                           %20.14e\n", END_T            );
      fprintf( Note, "END_STEP                        %ld\n",     END_STEP         );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of test problems
      fprintf( Note, "Parameters of Test Problems\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "TESTPROB_ID                     %d\n", TESTPROB_ID );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of code units
      fprintf( Note, "Parameters of Code Units\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__UNIT                       %d\n",                     OPT__UNIT                     );
      if ( OPT__UNIT ) {
#     ifdef COMOVING
      fprintf( Note, "\n### All units marked with (*) assume h = %13.7e ###\n\n", HUBBLE0 );

      const double Current_Matter_Density = OMEGA_M0*3*SQR( 100.0*HUBBLE0*Const_km/Const_Mpc/Const_s )/( 8.0*M_PI*Const_NewtonG );
      fprintf( Note, "rho_bg = current matter density = %20.14e Msun/kpc^3 (*)\n\n",
               Current_Matter_Density/(Const_Msun/CUBE(Const_kpc)) );

      fprintf( Note, "UNIT_L (length)                 %20.14e Mpc/h\n",          UNIT_L/(Const_Mpc/HUBBLE0)    );
      fprintf( Note, "                              = %20.14e cm         (*)\n", UNIT_L                        );
      fprintf( Note, "UNIT_M (mass)                   %20.14e Msun/h\n",         UNIT_M/(Const_Msun/HUBBLE0)   );
      fprintf( Note, "                              = %20.14e g          (*)\n", UNIT_M                        );
      fprintf( Note, "UNIT_T (time)                   %20.14e Gyr        (*)\n", UNIT_T/Const_Gyr              );
      fprintf( Note, "                              = %20.14e s          (*)\n", UNIT_T                        );
      fprintf( Note, "UNIT_V (velocity)               %20.14e km/s\n",           UNIT_V/(Const_km/Const_s)     );
      fprintf( Note, "                              = %20.14e cm/s\n",           UNIT_V                        );
      fprintf( Note, "UNIT_D (mass density)           %20.14e rho_bg     (*)\n", UNIT_D/Current_Matter_Density );
      fprintf( Note, "                              = %20.14e g/cm^3     (*)\n", UNIT_D                        );
      fprintf( Note, "UNIT_E (energy)                 %20.14e g*cm^2/s^2 (*)\n", UNIT_E                        );
      fprintf( Note, "UNIT_P (energy density)         %20.14e g/cm/s^2   (*)\n", UNIT_P                        );
#     ifdef MHD
#     error : ERROR : MHD is not supported here !!!
#     endif

#     else

      fprintf( Note, "UNIT_L                          %20.14e cm\n",             UNIT_L                        );
      fprintf( Note, "UNIT_M                          %20.14e g\n",              UNIT_M                        );
      fprintf( Note, "UNIT_T                          %20.14e s\n",              UNIT_T                        );
      fprintf( Note, "UNIT_V                          %20.14e cm/s\n",           UNIT_V                        );
      fprintf( Note, "UNIT_D                          %20.14e g/cm^3\n",         UNIT_D                        );
      fprintf( Note, "UNIT_E (energy)                 %20.14e g*cm^2/s^2\n",     UNIT_E                        );
      fprintf( Note, "UNIT_P (energy density)         %20.14e g/cm/s^2\n",       UNIT_P                        );
#     ifdef MHD
      fprintf( Note, "UNIT_B (magnetic field)         %20.14e gauss\n",          UNIT_B                        );
#     endif
#     endif // #ifdef COMOVING ... else ...
      }

      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of boundary condition
      fprintf( Note, "Parameters of Boundary Condition\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__BC_FLU[0] (-x)             %d\n",      OPT__BC_FLU[0] );
      fprintf( Note, "OPT__BC_FLU[1] (+x)             %d\n",      OPT__BC_FLU[1] );
      fprintf( Note, "OPT__BC_FLU[2] (-y)             %d\n",      OPT__BC_FLU[2] );
      fprintf( Note, "OPT__BC_FLU[3] (+y)             %d\n",      OPT__BC_FLU[3] );
      fprintf( Note, "OPT__BC_FLU[4] (-z)             %d\n",      OPT__BC_FLU[4] );
      fprintf( Note, "OPT__BC_FLU[5] (+z)             %d\n",      OPT__BC_FLU[5] );
#     ifdef GRAVITY
      fprintf( Note, "OPT__BC_POT                     %d\n",      OPT__BC_POT    );
      fprintf( Note, "GFUNC_COEFF0                    %13.7e\n",  GFUNC_COEFF0   );
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of particle
#     ifdef PARTICLE
      fprintf( Note, "Parameters of Particle\n" );
      fprintf( Note, "***********************************************************************************\n" );
#     ifdef DEBUG_PARTICLE
      fprintf( Note, "DEBUG_PARTICLE                  ON\n" );
#     else
      fprintf( Note, "DEBUG_PARTICLE                  OFF\n" );
#     endif
      fprintf( Note, "Par->NPar_Active_AllRank        %ld\n",     amr->Par->NPar_Active_AllRank );
      fprintf( Note, "Par->Init                       %d\n",      amr->Par->Init                );
      fprintf( Note, "Par->ParICFormat                %d\n",      amr->Par->ParICFormat         );
      fprintf( Note, "Par->ParICMass                 %14.7e\n",   amr->Par->ParICMass           );
      fprintf( Note, "Par->Interp                     %d\n",      amr->Par->Interp              );
      fprintf( Note, "Par->Integ                      %d\n",      amr->Par->Integ               );
      fprintf( Note, "Par->GhostSize                  %d\n",      amr->Par->GhostSize           );
      fprintf( Note, "Par->ImproveAcc                 %d\n",      amr->Par->ImproveAcc          );
      fprintf( Note, "Par->PredictPos                 %d\n",      amr->Par->PredictPos          );
      fprintf( Note, "Par->RemoveCell                 %13.7e\n",  amr->Par->RemoveCell          );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of cosmological simulations (comoving frame)
#     ifdef COMOVING
      fprintf( Note, "Parameters of Cosmological Simulation\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "A_INIT                          %13.7e\n",  A_INIT   );
      fprintf( Note, "OMEGA_M0                        %13.7e\n",  OMEGA_M0 );
      fprintf( Note, "HUBBLE0 (h)                     %13.7e\n",  HUBBLE0  );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of time-step determination
      fprintf( Note, "Parameters of Time-step Determination\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "DT__MAX                        %14.7e\n",   DT__MAX                   );
      fprintf( Note, "DT__FLUID                       %13.7e\n",  DT__FLUID                 );
      fprintf( Note, "DT__FLUID_INIT                  %13.7e\n",  DT__FLUID_INIT            );
#     ifdef SRHD
      fprintf( Note, "DT_SPEED_OF_LIGHT               %d\n",      DT_SPEED_OF_LIGHT         );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "DT__GRAVITY                     %13.7e\n",  DT__GRAVITY               );
#     endif
#     if ( MODEL == ELBDM )
      fprintf( Note, "DT__PHASE                       %13.7e\n",  DT__PHASE                 );
#     endif
#     ifdef PARTICLE
      fprintf( Note, "DT__PARVEL                      %13.7e\n",  DT__PARVEL                );
      fprintf( Note, "DT__PARVEL_MAX                 %14.7e\n",   DT__PARVEL_MAX            );
      fprintf( Note, "DT__PARACC                      %13.7e\n",  DT__PARACC                );
#     endif
#     ifdef COMOVING
      fprintf( Note, "DT__MAX_DELTA_A                 %13.7e\n",  DT__MAX_DELTA_A           );
#     endif
      fprintf( Note, "DT__SYNC_PARENT_LV              %13.7e\n",  DT__SYNC_PARENT_LV        );
      fprintf( Note, "DT__SYNC_CHILDREN_LV            %13.7e\n",  DT__SYNC_CHILDREN_LV      );
      fprintf( Note, "OPT__DT_USER                    %d\n",      OPT__DT_USER              );
      fprintf( Note, "OPT__DT_LEVEL                   %d\n",      OPT__DT_LEVEL             );
      fprintf( Note, "AUTO_REDUCE_DT                  %d\n",      AUTO_REDUCE_DT            );
      fprintf( Note, "AUTO_REDUCE_DT_FACTOR           %13.7e\n",  AUTO_REDUCE_DT_FACTOR     );
      fprintf( Note, "AUTO_REDUCE_DT_FACTOR_MIN       %13.7e\n",  AUTO_REDUCE_DT_FACTOR_MIN );
      fprintf( Note, "OPT__RECORD_DT                  %d\n",      OPT__RECORD_DT            );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of domain refinement
      fprintf( Note, "Parameters of Domain Refinement\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "REGRID_COUNT                    %d\n",      REGRID_COUNT              );
      fprintf( Note, "FLAG_BUFFER_SIZE                %d\n",      FLAG_BUFFER_SIZE          );
      fprintf( Note, "FLAG_BUFFER_SIZE_MAXM1_LV       %d\n",      FLAG_BUFFER_SIZE_MAXM1_LV );
      fprintf( Note, "FLAG_BUFFER_SIZE_MAXM2_LV       %d\n",      FLAG_BUFFER_SIZE_MAXM2_LV );
      fprintf( Note, "MAX_LEVEL                       %d\n",      MAX_LEVEL                 );
      fprintf( Note, "OPT__FLAG_RHO                   %d\n",      OPT__FLAG_RHO             );
      fprintf( Note, "OPT__FLAG_RHO_GRADIENT          %d\n",      OPT__FLAG_RHO_GRADIENT    );
#     if ( MODEL == HYDRO )
      fprintf( Note, "OPT__FLAG_PRES_GRADIENT         %d\n",      OPT__FLAG_PRES_GRADIENT   );
      fprintf( Note, "OPT__FLAG_ENGY_GRADIENT         %d\n",      OPT__FLAG_ENGY_GRADIENT   );
      fprintf( Note, "OPT__FLAG_VORTICITY             %d\n",      OPT__FLAG_VORTICITY       );
      fprintf( Note, "OPT__FLAG_JEANS                 %d\n",      OPT__FLAG_JEANS           );
#     ifdef MHD
      fprintf( Note, "OPT__FLAG_CURRENT               %d\n",      OPT__FLAG_CURRENT         );
#     endif
#     ifdef SRHD
      fprintf( Note, "OPT__FLAG_4VELOCITY             %d\n",      OPT__FLAG_4VELOCITY       );
      fprintf( Note, "OPT__FLAG_LORENTZ_GRADIENT      %d\n",      OPT__FLAG_LORENTZ_GRADIENT);
#     endif
#     endif
#     if ( MODEL == ELBDM )
      fprintf( Note, "OPT__FLAG_ENGY_DENSITY          %d\n",      OPT__FLAG_ENGY_DENSITY    );
#     endif
      fprintf( Note, "OPT__FLAG_LOHNER_DENS           %d\n",      OPT__FLAG_LOHNER_DENS     );
#     if ( MODEL == HYDRO )
      fprintf( Note, "OPT__FLAG_LOHNER_ENGY           %d\n",      OPT__FLAG_LOHNER_ENGY     );
      fprintf( Note, "OPT__FLAG_LOHNER_PRES           %d\n",      OPT__FLAG_LOHNER_PRES     );
      fprintf( Note, "OPT__FLAG_LOHNER_TEMP           %d\n",      OPT__FLAG_LOHNER_TEMP     );
#     endif
      fprintf( Note, "OPT__FLAG_LOHNER_FORM           %s\n",      (OPT__FLAG_LOHNER_FORM==LOHNER_FLASH1   ) ? "LOHNER_FLASH1"    :
                                                                  (OPT__FLAG_LOHNER_FORM==LOHNER_FLASH2   ) ? "LOHNER_FLASH2"    :
                                                                  (OPT__FLAG_LOHNER_FORM==LOHNER_FORM_INV1) ? "LOHNER_FORM_INV1" :
                                                                  (OPT__FLAG_LOHNER_FORM==LOHNER_FORM_INV2) ? "LOHNER_FORM_INV2" :
                                                                                                               "UNKNOWN" );
      fprintf( Note, "OPT__FLAG_USER                  %d\n",      OPT__FLAG_USER            );
      fprintf( Note, "OPT__FLAG_USER_NUM              %d\n",      OPT__FLAG_USER_NUM        );
      fprintf( Note, "OPT__FLAG_REGION                %d\n",      OPT__FLAG_REGION          );
#     ifdef PARTICLE
      fprintf( Note, "OPT__FLAG_NPAR_PATCH            %d\n",      OPT__FLAG_NPAR_PATCH      );
      fprintf( Note, "OPT__FLAG_NPAR_CELL             %d\n",      OPT__FLAG_NPAR_CELL       );
      fprintf( Note, "OPT__FLAG_PAR_MASS_CELL         %d\n",      OPT__FLAG_PAR_MASS_CELL   );
#     endif
      fprintf( Note, "OPT__NO_FLAG_NEAR_BOUNDARY      %d\n",      OPT__NO_FLAG_NEAR_BOUNDARY);
      fprintf( Note, "OPT__PATCH_COUNT                %d\n",      OPT__PATCH_COUNT          );
#     ifdef PARTICLE
      fprintf( Note, "OPT__PARTICLE_COUNT             %d\n",      OPT__PARTICLE_COUNT       );
#     endif
      fprintf( Note, "OPT__REUSE_MEMORY               %d\n",      OPT__REUSE_MEMORY         );
      fprintf( Note, "OPT__MEMORY_POOL                %d\n",      OPT__MEMORY_POOL          );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of parallelization
#     ifndef SERIAL
      fprintf( Note, "Parameters of Parallelization\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "Flu_ParaBuf                     %d\n",      Flu_ParaBuf               );
#     ifdef GRAVITY
      fprintf( Note, "Pot_ParaBuf                     %d\n",      Pot_ParaBuf               );
      fprintf( Note, "Rho_ParaBuf                     %d\n",      Rho_ParaBuf               );
#     endif
#     ifdef LOAD_BALANCE
      fprintf( Note, "LB_WLI_MAX                      %13.7e\n",  amr->LB->WLI_Max          );
#     ifdef PARTICLE
      fprintf( Note, "LB_PAR_WEIGHT                   %13.7e\n",  amr->LB->Par_Weight       );
#     endif
      fprintf( Note, "OPT__RECORD_LOAD_BALANCE        %d\n",      OPT__RECORD_LOAD_BALANCE  );
#     endif // #ifdef LOAD_BALANCE
      fprintf( Note, "OPT__MINIMIZE_MPI_BARRIER       %d\n",      OPT__MINIMIZE_MPI_BARRIER );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif // #ifndef SERIAL


//    record the parameters of source terms
      fprintf( Note, "Parameters of Source Terms\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "SRC_ANY                         %d\n",      SrcTerms.Any              );
      fprintf( Note, "SRC_DELEPTONIZATION             %d\n",      SrcTerms.Deleptonization  );
      fprintf( Note, "SRC_USER                        %d\n",      SrcTerms.User             );
      fprintf( Note, "SRC_GPU_NPGROUP                 %d\n",      SRC_GPU_NPGROUP           );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of Grackle
#     ifdef SUPPORT_GRACKLE
      fprintf( Note, "Parameters of Grackle\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "GRACKLE_ACTIVATE                %d\n",      GRACKLE_ACTIVATE        );
      if ( GRACKLE_ACTIVATE ) {
      fprintf( Note, "GRACKLE_VERBOSE                 %d\n",      GRACKLE_VERBOSE         );
      fprintf( Note, "GRACKLE_COOLING                 %d\n",      GRACKLE_COOLING         );
      fprintf( Note, "GRACKLE_PRIMORDIAL              %d\n",      GRACKLE_PRIMORDIAL      );
      fprintf( Note, "GRACKLE_METAL                   %d\n",      GRACKLE_METAL           );
      fprintf( Note, "GRACKLE_UV                      %d\n",      GRACKLE_UV              );
      fprintf( Note, "GRACKLE_CMB_FLOOR               %d\n",      GRACKLE_CMB_FLOOR       );
      fprintf( Note, "GRACKLE_PE_HEATING              %d\n",      GRACKLE_PE_HEATING      );
      fprintf( Note, "GRACKLE_PE_HEATING_RATE         %13.7e\n",  GRACKLE_PE_HEATING_RATE );
      fprintf( Note, "GRACKLE_CLOUDY_TABLE            %s\n",      GRACKLE_CLOUDY_TABLE    );
      fprintf( Note, "GRACKLE_THREE_BODY_RATE         %d\n",      GRACKLE_THREE_BODY_RATE );
      fprintf( Note, "GRACKLE_CIE_COOLING             %d\n",      GRACKLE_CIE_COOLING     );
      fprintf( Note, "GRACKLE_H2_OPA_APPROX           %d\n",      GRACKLE_H2_OPA_APPROX   );
      fprintf( Note, "CHE_GPU_NPGROUP                 %d\n",      CHE_GPU_NPGROUP         ); }
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif // #ifdef SUPPORT_GRACKLE


//    record the parameters of star formation
#     ifdef STAR_FORMATION
      fprintf( Note, "Parameters of Star Formation\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "SF_CREATE_STAR_SCHEME           %d\n",           SF_CREATE_STAR_SCHEME                          );
      if ( SF_CREATE_STAR_SCHEME != SF_CREATE_STAR_SCHEME_NONE ) {
      fprintf( Note, "SF_CREATE_STAR_RSEED            %d\n",           SF_CREATE_STAR_RSEED                           );
      fprintf( Note, "SF_CREATE_STAR_DET_RANDOM       %d\n",           SF_CREATE_STAR_DET_RANDOM                      );
      fprintf( Note, "SF_CREATE_STAR_MIN_LEVEL        %d\n",           SF_CREATE_STAR_MIN_LEVEL                       );
      fprintf( Note, "SF_CREATE_STAR_MIN_GAS_DENS     %13.7e\n",       SF_CREATE_STAR_MIN_GAS_DENS                    );
      fprintf( Note, "                              = %13.7e cm^-3\n", SF_CREATE_STAR_MIN_GAS_DENS*UNIT_D/Const_mH    );
      fprintf( Note, "SF_CREATE_STAR_MASS_EFF         %13.7e\n",       SF_CREATE_STAR_MASS_EFF                        );
      fprintf( Note, "SF_CREATE_STAR_MIN_STAR_MASS    %13.7e\n",       SF_CREATE_STAR_MIN_STAR_MASS                   );
      fprintf( Note, "                              = %13.7e Msun\n",  SF_CREATE_STAR_MIN_STAR_MASS*UNIT_M/Const_Msun );
      fprintf( Note, "SF_CREATE_STAR_MAX_STAR_MFRAC   %13.7e\n",       SF_CREATE_STAR_MAX_STAR_MFRAC                  ); }
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif // #ifdef STAR_FORMATION


//    record the parameters of Fluid solver in different models
      fprintf( Note, "Parameters of Fluid Solver (in different models)\n" );
      fprintf( Note, "***********************************************************************************\n" );
#     if   ( MODEL == HYDRO )
      fprintf( Note, "GAMMA                           %13.7e\n",  GAMMA                   );
      fprintf( Note, "MOLECULAR_WEIGHT                %13.7e\n",  MOLECULAR_WEIGHT        );
      fprintf( Note, "ISO_TEMP                        %13.7e\n",  ISO_TEMP                );
      fprintf( Note, "MINMOD_COEFF                    %13.7e\n",  MINMOD_COEFF            );
      fprintf( Note, "OPT__LR_LIMITER                 %s\n",      ( OPT__LR_LIMITER == VANLEER           ) ? "VANLEER"    :
                                                                  ( OPT__LR_LIMITER == GMINMOD           ) ? "GMINMOD"    :
                                                                  ( OPT__LR_LIMITER == ALBADA            ) ? "ALBADA"     :
                                                                  ( OPT__LR_LIMITER == VL_GMINMOD        ) ? "VL_GMINMOD" :
                                                                  ( OPT__LR_LIMITER == EXTPRE            ) ? "EXTPRE"     :
                                                                  ( OPT__LR_LIMITER == LR_LIMITER_NONE   ) ? "NONE"       :
                                                                                                             "UNKNOWN" );
      fprintf( Note, "OPT__1ST_FLUX_CORR              %s\n",      ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D   ) ? "3D"   :
                                                                  ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_3D1D ) ? "3D1D" :
                                                                  ( OPT__1ST_FLUX_CORR == FIRST_FLUX_CORR_NONE ) ? "NONE" :
                                                                                                                   "UNKNOWN" );
      fprintf( Note, "OPT__1ST_FLUX_CORR_SCHEME       %s\n",      ( OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_ROE  ) ? "RSOLVER_1ST_ROE"  :
                                                                  ( OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_HLLC ) ? "RSOLVER_1ST_HLLC" :
                                                                  ( OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_HLLE ) ? "RSOLVER_1ST_HLLE" :
                                                                  ( OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_HLLD ) ? "RSOLVER_1ST_HLLD" :
                                                                  ( OPT__1ST_FLUX_CORR_SCHEME == RSOLVER_1ST_NONE ) ? "NONE"             :
                                                                                                                "UNKNOWN" );

#     elif ( MODEL == ELBDM )
      if ( OPT__UNIT ) {
//    since the mass unit in cosmological simulation has the 1/h dependence, the actual ELBDM_MASS adopted in the
//    cosmological simulation also depends on 1/h
//    --> however, Planck constant also depends on 1/h^2 --> ELBDM_ETA depends on h
//    --> since Planck constant is a real constant which cannot be changed, it's more appropriate to express
//        ELBDM particle mass as ev/c^2*h (not ev/c^2/h), just like the length unit Mpc/h
//    --> for a given simulation result, we can always reinterpret h to give different box size and ELBDM particle mass
//    --> for example, the simulations results with h=1, m=1.0e-22 eV can be reinterpreted as h=0.5, m=5.0e-23 eV
//    --> also note that this data reinterpretation is purely based on redefining basic units and is different from the
//        scaling symmetry in ELBDM
#     ifdef COMOVING
      fprintf( Note, "ELBDM_MASS                      %13.7e %s\n",  ELBDM_MASS*UNIT_M/(HUBBLE0*Const_eV/SQR(Const_c)), "h*ev/c^2" );
      fprintf( Note, "                              = %13.7e %s (assuming h = %13.7e)\n",
                                                                     ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c)), "ev/c^2", HUBBLE0 );
#     else
      fprintf( Note, "ELBDM_MASS                      %13.7e %s\n",  ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c)), "ev/c^2" );
#     endif
      }
      else
      fprintf( Note, "ELBDM_MASS                      %13.7e\n",     ELBDM_MASS           );
      fprintf( Note, "ELBDM_PLANCK_CONST              %13.7e\n",     ELBDM_PLANCK_CONST   );
      fprintf( Note, "ELBDM_ETA                       %13.7e\n",     ELBDM_ETA            );
#     ifdef QUARTIC_SELF_INTERACTION
      fprintf( Note, "ELBDM_LAMBDA                    %13.7e\n",     ELBDM_LAMBDA         );
#     endif
      fprintf( Note, "ELBDM_TAYLOR3_COEFF             %13.7e\n",     ELBDM_TAYLOR3_COEFF  );
      fprintf( Note, "ELBDM_TAYLOR3_AUTO              %d\n",         ELBDM_TAYLOR3_AUTO   );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of Fluid solver in different models
      fprintf( Note, "Parameters of Fluid Solver (in all models)\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "FLU_GPU_NPGROUP                 %d\n",      FLU_GPU_NPGROUP          );
      fprintf( Note, "GPU_NSTREAM                     %d\n",      GPU_NSTREAM              );
      fprintf( Note, "OPT__FIXUP_FLUX                 %d\n",      OPT__FIXUP_FLUX          );
#     ifdef MHD
      fprintf( Note, "OPT__FIXUP_ELECTRIC             %d\n",      OPT__FIXUP_ELECTRIC      );
#     endif
      fprintf( Note, "OPT__FIXUP_RESTRICT             %d\n",      OPT__FIXUP_RESTRICT      );
      fprintf( Note, "OPT__CORR_AFTER_ALL_SYNC        %d\n",      OPT__CORR_AFTER_ALL_SYNC );
      fprintf( Note, "OPT__NORMALIZE_PASSIVE          %d\n",      OPT__NORMALIZE_PASSIVE   );

//    target passive scalars to be normalized
      fprintf( Note, "   Number of fields             %d\n",      PassiveNorm_NVar         );
      if ( PassiveNorm_NVar > 0 ) {
      fprintf( Note, "   Target fields               "                                     );
      for (int v=0; v<PassiveNorm_NVar; v++)
      fprintf( Note, " %s",                                       FieldLabel[ NCOMP_FLUID + PassiveNorm_VarIdx[v] ] );
      fprintf( Note, "\n" ); }

      fprintf( Note, "OPT__OVERLAP_MPI                %d\n",      OPT__OVERLAP_MPI         );
      fprintf( Note, "OPT__RESET_FLUID                %d\n",      OPT__RESET_FLUID         );
#     if ( MODEL == HYDRO  ||  MODEL == ELBDM )
      fprintf( Note, "MIN_DENS                        %13.7e\n",  MIN_DENS                 );
#     endif
#     if ( MODEL == HYDRO )
      fprintf( Note, "MIN_PRES                        %13.7e\n",  MIN_PRES                 );
      fprintf( Note, "MIN_EINT                        %13.7e\n",  MIN_EINT                 );
      fprintf( Note, "OPT__LAST_RESORT_FLOOR          %d\n",      OPT__LAST_RESORT_FLOOR   );
      fprintf( Note, "JEANS_MIN_PRES                  %d\n",      JEANS_MIN_PRES           );
      if ( JEANS_MIN_PRES ) {
      fprintf( Note, "JEANS_MIN_PRES_LEVEL            %d\n",      JEANS_MIN_PRES_LEVEL     );
      fprintf( Note, "JEANS_MIN_PRES_NCELL            %d\n",      JEANS_MIN_PRES_NCELL     ); }
#     endif
#     ifdef DUAL_ENERGY
      fprintf( Note, "DUAL_ENERGY_SWITCH              %13.7e\n",  DUAL_ENERGY_SWITCH       );
#     endif
      fprintf( Note, "WITH_COARSE_FINE_FLUX           %d\n",      amr->WithFlux            );
#     ifdef MHD
      fprintf( Note, "WITH_COARSE_FINE_ELECTRIC       %d\n",      amr->WithElectric        );
#     endif
#     ifndef SERIAL
      int MPI_Thread_Status;
      MPI_Query_thread( &MPI_Thread_Status );
      fprintf( Note, "MPI Thread Level                " );
      switch ( MPI_Thread_Status )
      {
         case MPI_THREAD_SINGLE:       fprintf( Note, "MPI_THREAD_SINGLE\n" );       break;
         case MPI_THREAD_FUNNELED:     fprintf( Note, "MPI_THREAD_FUNNELED\n" );     break;
         case MPI_THREAD_SERIALIZED:   fprintf( Note, "MPI_THREAD_SERIALIZED\n" );   break;
         case MPI_THREAD_MULTIPLE:     fprintf( Note, "MPI_THREAD_MULTIPLE\n" );     break;

         default:                      fprintf( Note, "UNKNOWN\n" );
      }
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of Poisson and Gravity solvers
#     ifdef GRAVITY
      fprintf( Note, "Parameters of Poisson and Gravity Solvers\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "NEWTON_G                        %13.7e\n",  NEWTON_G                );
#     if   ( POT_SCHEME == SOR )
      fprintf( Note, "SOR_OMEGA                       %13.7e\n",  SOR_OMEGA               );
      fprintf( Note, "SOR_MAX_ITER                    %d\n",      SOR_MAX_ITER            );
      fprintf( Note, "SOR_MIN_ITER                    %d\n",      SOR_MIN_ITER            );
#     elif ( POT_SCHEME == MG )
      fprintf( Note, "MG_MAX_ITER                     %d\n",      MG_MAX_ITER             );
      fprintf( Note, "MG_NPRE_SMOOTH                  %d\n",      MG_NPRE_SMOOTH          );
      fprintf( Note, "MG_NPOST_SMOOTH                 %d\n",      MG_NPOST_SMOOTH         );
      fprintf( Note, "MG_TOLERATED_ERROR              %13.7e\n",  MG_TOLERATED_ERROR      );
#     endif
      fprintf( Note, "POT_GPU_NPGROUP                 %d\n",      POT_GPU_NPGROUP         );
      fprintf( Note, "OPT__GRA_P5_GRADIENT            %d\n",      OPT__GRA_P5_GRADIENT    );
      fprintf( Note, "OPT__SELF_GRAVITY               %d\n",      OPT__SELF_GRAVITY       );
      fprintf( Note, "OPT__EXT_ACC                    %d\n",      OPT__EXT_ACC            );
      fprintf( Note, "OPT__EXT_POT                    %d\n",      OPT__EXT_POT            );
      if ( OPT__EXT_POT == EXT_POT_TABLE ) {
      fprintf( Note, "EXT_POT_TABLE_NAME              %s\n",      EXT_POT_TABLE_NAME      );
      fprintf( Note, "EXT_POT_TABLE_NPOINT_X          %d\n",      EXT_POT_TABLE_NPOINT[0] );
      fprintf( Note, "EXT_POT_TABLE_NPOINT_Y          %d\n",      EXT_POT_TABLE_NPOINT[1] );
      fprintf( Note, "EXT_POT_TABLE_NPOINT_Z          %d\n",      EXT_POT_TABLE_NPOINT[2] );
      fprintf( Note, "EXT_POT_TABLE_DH                %13.7e\n",  EXT_POT_TABLE_DH        );
      fprintf( Note, "EXT_POT_TABLE_EDGEL_X          %14.7e\n",   EXT_POT_TABLE_EDGEL[0]  );
      fprintf( Note, "EXT_POT_TABLE_EDGEL_Y          %14.7e\n",   EXT_POT_TABLE_EDGEL[1]  );
      fprintf( Note, "EXT_POT_TABLE_EDGEL_Z          %14.7e\n",   EXT_POT_TABLE_EDGEL[2]  );
      fprintf( Note, "EXT_POT_TABLE_FLOAT8            %d\n",      EXT_POT_TABLE_FLOAT8    ); }
      fprintf( Note, "OPT__GRAVITY_EXTRA_MASS         %d\n",      OPT__GRAVITY_EXTRA_MASS );
      fprintf( Note, "AveDensity_Init                 %13.7e\n",  AveDensity_Init         );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif // #ifdef GRAVITY


//    record the parameters of initialization
      fprintf( Note, "Parameters of Initialization\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__INIT                       %d\n",      OPT__INIT               );
      fprintf( Note, "RESTART_LOAD_NRANK              %d\n",      RESTART_LOAD_NRANK      );
      fprintf( Note, "OPT__RESTART_RESET              %d\n",      OPT__RESTART_RESET      );
      fprintf( Note, "OPT__UM_IC_LEVEL                %d\n",      OPT__UM_IC_LEVEL        );
      fprintf( Note, "OPT__UM_IC_NVAR                 %d\n",      OPT__UM_IC_NVAR         );
      fprintf( Note, "OPT__UM_IC_FORMAT               %d\n",      OPT__UM_IC_FORMAT       );
      fprintf( Note, "OPT__UM_IC_DOWNGRADE            %d\n",      OPT__UM_IC_DOWNGRADE    );
      fprintf( Note, "OPT__UM_IC_REFINE               %d\n",      OPT__UM_IC_REFINE       );
      fprintf( Note, "OPT__UM_IC_LOAD_NRANK           %d\n",      OPT__UM_IC_LOAD_NRANK   );
      fprintf( Note, "OPT__INIT_RESTRICT              %d\n",      OPT__INIT_RESTRICT      );
      fprintf( Note, "OPT__INIT_GRID_WITH_OMP         %d\n",      OPT__INIT_GRID_WITH_OMP );
      fprintf( Note, "OPT__GPUID_SELECT               %d\n",      OPT__GPUID_SELECT       );
      fprintf( Note, "INIT_SUBSAMPLING_NCELL          %d\n",      INIT_SUBSAMPLING_NCELL  );
#     ifdef MHD
      fprintf( Note, "OPT__INIT_BFIELD_BYFILE         %d\n",      OPT__INIT_BFIELD_BYFILE );
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of interpolation schemes
      fprintf( Note, "Parameters of Interpolation Schemes\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__INT_TIME                   %d\n",      OPT__INT_TIME           );
#     if ( MODEL == ELBDM )
      fprintf( Note, "OPT__INT_PHASE                  %d\n",      OPT__INT_PHASE          );
#     endif
      fprintf( Note, "OPT__FLU_INT_SCHEME             %s\n",      ( OPT__FLU_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                                  ( OPT__FLU_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                                  ( OPT__FLU_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                                  ( OPT__FLU_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                                  ( OPT__FLU_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                                  ( OPT__FLU_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                                  ( OPT__FLU_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                            "UNKNOWN" );
#     ifdef MHD
      fprintf( Note, "OPT__MAG_INT_SCHEME             %s\n",      ( OPT__MAG_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                                  ( OPT__MAG_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                                  ( OPT__MAG_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                                  ( OPT__MAG_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                                  ( OPT__MAG_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                                  ( OPT__MAG_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                                  ( OPT__MAG_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                            "UNKNOWN" );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "OPT__POT_INT_SCHEME             %s\n",      ( OPT__POT_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                                  ( OPT__POT_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                                  ( OPT__POT_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                                  ( OPT__POT_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                                  ( OPT__POT_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                                  ( OPT__POT_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                                  ( OPT__POT_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                            "UNKNOWN" );
      fprintf( Note, "OPT__RHO_INT_SCHEME             %s\n",      ( OPT__RHO_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                                  ( OPT__RHO_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                                  ( OPT__RHO_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                                  ( OPT__RHO_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                                  ( OPT__RHO_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                                  ( OPT__RHO_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                                  ( OPT__RHO_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                            "UNKNOWN" );
      fprintf( Note, "OPT__GRA_INT_SCHEME             %s\n",      ( OPT__GRA_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                                  ( OPT__GRA_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                                  ( OPT__GRA_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                                  ( OPT__GRA_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                                  ( OPT__GRA_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                                  ( OPT__GRA_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                                  ( OPT__GRA_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                            "UNKNOWN" );
#     endif
      fprintf( Note, "OPT__REF_FLU_INT_SCHEME         %s\n",   ( OPT__REF_FLU_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                               ( OPT__REF_FLU_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                             "UNKNOWN" );
#     ifdef MHD
      fprintf( Note, "OPT__REF_MAG_INT_SCHEME         %s\n",   ( OPT__REF_MAG_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                               ( OPT__REF_MAG_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                             "UNKNOWN" );
#     endif
#     ifdef GRAVITY
      fprintf( Note, "OPT__REF_POT_INT_SCHEME         %s\n",   ( OPT__REF_POT_INT_SCHEME == INT_MINMOD3D ) ? "MINMOD3D" :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_MINMOD1D ) ? "MINMOD1D" :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_VANLEER  ) ? "VANLEER"  :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_CQUAD    ) ? "CQUAD"    :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_QUAD     ) ? "QUAD"     :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_CQUAR    ) ? "CQUAR"    :
                                                               ( OPT__REF_POT_INT_SCHEME == INT_QUAR     ) ? "QUAR"     :
                                                                                                             "UNKNOWN" );
#     endif
      fprintf( Note, "INT_MONO_COEFF                  %13.7e\n",  INT_MONO_COEFF          );
      fprintf( Note, "INT_OPP_SIGN_0TH_ORDER          %d\n",      INT_OPP_SIGN_0TH_ORDER  );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of data dump
      fprintf( Note, "Parameters of Data Dump\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__OUTPUT_TOTAL               %d\n",      OPT__OUTPUT_TOTAL    );
      fprintf( Note, "OPT__OUTPUT_PART                %d\n",      OPT__OUTPUT_PART     );
      fprintf( Note, "OPT__OUTPUT_USER                %d\n",      OPT__OUTPUT_USER     );
#     ifdef PARTICLE
      fprintf( Note, "OPT__OUTPUT_PAR_TEXT            %d\n",      OPT__OUTPUT_PAR_TEXT );
#     endif
      fprintf( Note, "OPT__OUTPUT_BASEPS              %d\n",      OPT__OUTPUT_BASEPS   );
      fprintf( Note, "OPT__OUTPUT_BASE                %d\n",      OPT__OUTPUT_BASE     );
#     ifdef GRAVITY
      fprintf( Note, "OPT__OUTPUT_POT                 %d\n",      OPT__OUTPUT_POT      );
#     endif
#     ifdef PARTICLE
      fprintf( Note, "OPT__OUTPUT_PAR_DENS            %d\n",      OPT__OUTPUT_PAR_DENS );
#     endif
#     ifdef MHD
      fprintf( Note, "OPT__OUTPUT_CC_MAG              %d\n",      OPT__OUTPUT_CC_MAG   );
#     endif
      fprintf( Note, "OPT__OUTPUT_MODE                %d\n",      OPT__OUTPUT_MODE     );
      fprintf( Note, "OUTPUT_STEP                     %d\n",      OUTPUT_STEP          );
      fprintf( Note, "OUTPUT_DT                       %20.14e\n", OUTPUT_DT            );
      fprintf( Note, "OUTPUT_PART_X                   %20.14e\n", OUTPUT_PART_X        );
      fprintf( Note, "OUTPUT_PART_Y                   %20.14e\n", OUTPUT_PART_Y        );
      fprintf( Note, "OUTPUT_PART_Z                   %20.14e\n", OUTPUT_PART_Z        );
      fprintf( Note, "INIT_DUMPID                     %d\n",      INIT_DUMPID          );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of yt inline analysis
#     ifdef SUPPORT_LIBYT
      fprintf( Note, "Parameters of YT Inline Analysis\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "YT_SCRIPT                       %s\n",      YT_SCRIPT  );
      fprintf( Note, "YT_VERBOSE                      %d\n",      YT_VERBOSE );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
#     endif


//    record the parameters of miscellaneous purposes
      fprintf( Note, "Parameters of Miscellaneous Purposes\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__VERBOSE                    %d\n",      OPT__VERBOSE             );
      fprintf( Note, "OPT__TIMING_BARRIER             %d\n",      OPT__TIMING_BARRIER      );
      fprintf( Note, "OPT__TIMING_BALANCE             %d\n",      OPT__TIMING_BALANCE      );
      fprintf( Note, "OPT__TIMING_MPI                 %d\n",      OPT__TIMING_MPI          );
      fprintf( Note, "OPT__RECORD_NOTE                %d\n",      OPT__RECORD_NOTE         );
      fprintf( Note, "OPT__RECORD_UNPHY               %d\n",      OPT__RECORD_UNPHY        );
      fprintf( Note, "OPT__RECORD_MEMORY              %d\n",      OPT__RECORD_MEMORY       );
      fprintf( Note, "OPT__RECORD_PERFORMANCE         %d\n",      OPT__RECORD_PERFORMANCE  );
      fprintf( Note, "OPT__MANUAL_CONTROL             %d\n",      OPT__MANUAL_CONTROL      );
      fprintf( Note, "OPT__RECORD_USER                %d\n",      OPT__RECORD_USER         );
      fprintf( Note, "OPT__OPTIMIZE_AGGRESSIVE        %d\n",      OPT__OPTIMIZE_AGGRESSIVE );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the parameters of simulation checks
      fprintf( Note, "Parameters of Simulation Checks\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OPT__CK_REFINE                  %d\n",      OPT__CK_REFINE            );
      fprintf( Note, "OPT__CK_PROPER_NESTING          %d\n",      OPT__CK_PROPER_NESTING    );
      fprintf( Note, "OPT__CK_CONSERVATION            %d\n",      OPT__CK_CONSERVATION      );
      fprintf( Note, "OPT__CK_NORMALIZE_PASSIVE       %d\n",      OPT__CK_NORMALIZE_PASSIVE );
      fprintf( Note, "OPT__CK_RESTRICT                %d\n",      OPT__CK_RESTRICT          );
      fprintf( Note, "OPT__CK_FINITE                  %d\n",      OPT__CK_FINITE            );
      fprintf( Note, "OPT__CK_PATCH_ALLOCATE          %d\n",      OPT__CK_PATCH_ALLOCATE    );
      fprintf( Note, "OPT__CK_FLUX_ALLOCATE           %d\n",      OPT__CK_FLUX_ALLOCATE     );
#     if ( MODEL == HYDRO )
      fprintf( Note, "OPT__CK_NEGATIVE                %d\n",      OPT__CK_NEGATIVE          );
#     endif
      fprintf( Note, "OPT__CK_MEMFREE                 %13.7e\n",  OPT__CK_MEMFREE           );
#     ifdef PARTICLE
      fprintf( Note, "OPT__CK_PARTICLE                %d\n",      OPT__CK_PARTICLE          );
#     endif
#     ifdef MHD
      fprintf( Note, "OPT__CK_INTERFACE_B             %d\n",      OPT__CK_INTERFACE_B       );
      fprintf( Note, "OPT__CK_DIVERGENCE_B            %d\n",      OPT__CK_DIVERGENCE_B      );
#     endif
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the flag criterion (density/density gradient/pressure gradient/user-defined)
      if ( OPT__FLAG_RHO )
      {
         fprintf( Note, "Flag Criterion (Density)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level             Density\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Rho[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_RHO_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Density Gradient)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level    Density Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_RhoGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_PRES_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Pressure Gradient in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Pressure Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_PresGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_ENGY_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Reduced Energy Density Gradient in SRHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Reduced Energy Density Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_EngyGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     ifdef SRHD
      if ( OPT__FLAG_LORENTZ_GRADIENT )
      {
         fprintf( Note, "Flag Criterion (Lorentz factor Gradient in SRHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Lorentz factor Gradient\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_LorentzFactorGradient[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_MOM_OVER_DENS )
      {
         fprintf( Note, "Flag Criterion ( Momentum over density in SRHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level   Momentum over density\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Mom_Over_Dens[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     else
      if ( OPT__FLAG_VORTICITY )
      {
         fprintf( Note, "Flag Criterion (Vorticity in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level           Vorticity\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Vorticity[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_JEANS )
      {
         fprintf( Note, "Flag Criterion (Jeans Length over Cell Size in HYDRO)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level       lambda_J / dh\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Jeans[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif

#     ifdef MHD
      if ( OPT__FLAG_CURRENT )
      {
         fprintf( Note, "Flag Criterion (Current Density in MHD)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level             Current\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_Current[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif
#     endif // #if ( MODEL == HYDRO )

#     if ( MODEL == ELBDM )
      if ( OPT__FLAG_ENGY_DENSITY )
      {
         fprintf( Note, "Flag Criterion (Energy Density in ELBDM)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level     Angle_over_2*PI              Soften\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
            fprintf( Note, "%7d%20.7e%20.7e\n", lv, FlagTable_EngyDensity[lv][0], FlagTable_EngyDensity[lv][1] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif

#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES || OPT__FLAG_LOHNER_TEMP )
#     elif ( MODEL == ELBDM )
      if ( OPT__FLAG_LOHNER_DENS )
#     endif
      {
         fprintf( Note, "Flag Criterion (Lohner Error Estimator)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level           Threshold              Filter              Soften      MinimumDensity\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
            fprintf( Note, "%7d%20.7e%20.7e%20.7e%20.7e\n", lv, FlagTable_Lohner[lv][0], FlagTable_Lohner[lv][1],
                     FlagTable_Lohner[lv][2], FlagTable_Lohner[lv][3] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_USER )
      {
         fprintf( Note, "Flag Criterion (User-defined)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level           Threshold\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)
         {
                                                       fprintf( Note, "%7d",    lv );
            for (int t=0; t<OPT__FLAG_USER_NUM; t++)   fprintf( Note, "%20.7e", FlagTable_User[lv][t] );
                                                       fprintf( Note, "\n" );
         }
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

#     ifdef PARTICLE
      if ( OPT__FLAG_NPAR_PATCH )
      {
         fprintf( Note, "Flag Criterion (# of Particles per Patch)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level      # of Particles\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20d\n", lv, FlagTable_NParPatch[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_NPAR_CELL )
      {
         fprintf( Note, "Flag Criterion (# of Particles per Cell)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level      # of Particles\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20d\n", lv, FlagTable_NParCell[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }

      if ( OPT__FLAG_PAR_MASS_CELL )
      {
         fprintf( Note, "Flag Criterion (Particle Mass per Cell)\n" );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "  Level       Particle Mass\n" );
         for (int lv=0; lv<MAX_LEVEL; lv++)  fprintf( Note, "%7d%20.7e\n", lv, FlagTable_ParMassCell[lv] );
         fprintf( Note, "***********************************************************************************\n" );
         fprintf( Note, "\n\n");
      }
#     endif // #ifdef PARTICLE


//    record the grid size in different refinement level
      fprintf( Note, "Cell Size and Scale (scale = number of cells at the finest level)\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "%7s%*c%26s%*c%16s\n", "Level", 10, ' ', "Cell Size", 10, ' ', "Cell Scale" );
      for (int lv=0; lv<NLEVEL; lv++)
      fprintf( Note, "%7d%*c%26.20lf%*c%16d\n", lv, 10, ' ', amr->dh[lv], 10, ' ', amr->scale[lv] );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


//    record the compilation time (of the file "Aux_TakeNote")
      fprintf( Note, "Compilation Time\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "%s %s\n", __DATE__, __TIME__ );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");


      fclose( Note );
   } // if ( MPI_Rank == 0 )


// record the hostname and PID of each MPI process (the function "CUAPI_DiagnoseDevice" will also record them)
#  ifndef GPU
   const int PID = getpid();
   char Host[1024];
   gethostname( Host, 1024 );

   if ( MPI_Rank == 0 )
   {
       Note = fopen( FileName, "a" );
       fprintf( Note, "Device Diagnosis\n" );
       fprintf( Note, "***********************************************************************************\n" );
       fclose( Note );
   }

   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         Note = fopen( FileName, "a" );
         if ( MPI_Rank != 0 )    fprintf( Note, "\n" );
         fprintf( Note, "MPI_Rank = %3d, hostname = %10s, PID = %5d\n", MPI_Rank, Host, PID );
         fprintf( Note, "CPU Info :\n" );
         fflush( Note );

         Aux_GetCPUInfo( FileName );

         fclose( Note );
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
      fclose( Note );
   }
#  endif // #ifndef GPU


// record the OpenMP status
#  ifdef OPENMP
   int omp_nthread, omp_chunk_size, omp_nested;
   omp_sched_t omp_schedule;

   omp_nested = omp_get_nested();
   omp_get_schedule( &omp_schedule, &omp_chunk_size );

#  pragma omp parallel
#  pragma omp master
   { omp_nthread = omp_get_num_threads(); }

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "OpenMP Diagnosis\n" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "OMP__SCHEDULE                   %s\n",      ( omp_schedule == omp_sched_static  ) ? "STATIC"  :
                                                                  ( omp_schedule == omp_sched_dynamic ) ? "DYNAMIC" :
                                                                  ( omp_schedule == omp_sched_guided  ) ? "GUIDED"  :
                                                                  ( omp_schedule == omp_sched_auto    ) ? "AUTO"    : "UNKNOWN" );
      fprintf( Note, "OMP__SCHEDULE_CHUNK_SIZE        %d\n",      omp_chunk_size          );
      fprintf( Note, "OMP__NESTED                     %s\n",      ( omp_nested ) ? "ON" : "OFF" );
      fprintf( Note, "\n" );
      fprintf( Note, "CPU core IDs of all OpenMP threads (tid == thread ID):\n" );
      fprintf( Note, "------------------------------------------------------------------------\n" );
      fprintf( Note, "%5s  %10s  %7s", "Rank", "Host", "NThread" );
      for (int t=0; t<omp_nthread; t++)
      fprintf( Note, "  tid-%02d", t );
      fprintf( Note, "\n" );
      fclose( Note );
   }

// record the CPU core id of each OpenMP thread
   int *omp_core_id = new int [omp_nthread];

#  pragma omp parallel
   { omp_core_id[ omp_get_thread_num() ] = get_cpuid(); }

   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         char MPI_Host[1024];
         gethostname( MPI_Host, 1024 );

         Note = fopen( FileName, "a" );
         fprintf( Note, "%5d  %10s  %7d", MPI_Rank, MPI_Host, omp_nthread );
         for (int t=0; t<omp_nthread; t++)   fprintf( Note, "  %6d", omp_core_id[t] );
         fprintf( Note, "\n" );
         fflush( Note );
         fclose( Note );
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }

   delete [] omp_core_id;

   if ( MPI_Rank == 0 )
   {
      Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fprintf( Note, "\n\n");
      fclose( Note );
   }
#  endif // #ifdef OPENMP


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_TakeNote ... done\n" );

} // FUNCTION : Aux_TakeNote



//-------------------------------------------------------------------------------------------------------
// Function    :  get_cpuid
// Description :  Get the CPU ID
//
// Note        :  Work on both macOS and Linux systems
//-------------------------------------------------------------------------------------------------------
int get_cpuid()
{

// See https://stackoverflow.com/questions/33745364/sched-getcpu-equivalent-for-os-x
   int CPU;

#  ifdef __APPLE__
   uint32_t CPUInfo[4];
   __cpuid_count(1, 0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
   if ((CPUInfo[3] & (1 << 9)) == 0) {
      CPU = -1;  /* no APIC on chip */
   } else {
      CPU = (unsigned)CPUInfo[1] >> 24;
   }
   if (CPU < 0) CPU = 0;
#  else
   CPU = sched_getcpu();
#  endif

   return CPU;

} // FUNCTION : get_cpuid
