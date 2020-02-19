#ifndef __MACRO_H__
#define __MACRO_H__



// ****************************************************************************
// ** This header defines the symbolic constants and macros used in GAMER.   **
// ** For clarity, useless options defined in the makefile will be "undef"   **
// ** in the end of this file.                                               **
// ****************************************************************************


// ########################
// ## Symbolic Constants ##
// ########################

// current version
#define VERSION      "gamer-2.1.1.dev"


// option == NONE --> the option is turned off
#define NONE         0


// GPU architecture
#define FERMI        1
#define KEPLER       2
#define MAXWELL      3
#define PASCAL       4
#define VOLTA        5
#define TURING       6


// models
#define HYDRO        1
//#define MHD        2     // MHD is now regarded as an option of HYDRO
#define ELBDM        3
#define PAR_ONLY     4


// hydrodynamic schemes
#define RTVD         1
#define MHM          3
#define MHM_RP       4
#define CTU          5


// data reconstruction schemes
#define PLM          1
#define PPM          2


// Riemann solvers
#define EXACT        1
#define ROE          2
#define HLLE         3
#define HLLC         4
#define HLLD         5


// dual-energy variables
#define DE_ENPY      1
#define DE_EINT      2

#ifdef DUAL_ENERGY
#define DE_UPDATED_BY_ETOT       ('0')
#define DE_UPDATED_BY_DUAL       ('1')
#define DE_UPDATED_BY_MIN_PRES   ('2')
#define DE_UPDATED_BY_ETOT_GRA   ('3')
#endif


// equation of states
#define EOS_GAMMA       1
#define EOS_ISOTHERMAL  2
#define EOS_NUCLEAR     3
#define EOS_TABULAR     4
#define EOS_USER        5


// Poisson solvers
#define SOR          1
#define MG           2


// load-balance parallelization
#define HILBERT      1


// random number implementation
#define RNG_GNU_EXT  1
#define RNG_CPP11    2


// NCOMP_FLUID : number of active components in each cell (for patch->fluid[])
//               --> do not include passive components here, which is set by NCOMP_PASSIVE
// NFLUX_FLUID : number of active components in patch->flux[]
//               --> do not include passive components here, which is set by NFLUX_PASSIVE
// NCOMP_MAG   : number of magnetic field components (for patch->magnetic[])
// NCOMP_ELE   : number of electric field components on each cell face (for patch->electric[])
#if   ( MODEL == HYDRO )
#  define NCOMP_FLUID         5
#  define NFLUX_FLUID         NCOMP_FLUID
# ifdef MHD
#  define NCOMP_MAG           3
#  define NCOMP_ELE           2
# else
#  define NCOMP_MAG           0
#  define NCOMP_ELE           0
# endif

// for ELBDM, we only need the density flux
#elif ( MODEL == ELBDM )
#  define NCOMP_FLUID         3
#  define NFLUX_FLUID         1

#elif ( MODEL == PAR_ONLY )
#  define NCOMP_FLUID         0
#  define NFLUX_FLUID         0

#else
#  error : ERROR : unsupported MODEL (please edit NCOMP_FLUID and NFLUX_FLUID for the new MODEL) !!
#endif // MODEL


// number of passively advected components in each cell

// define NCOMP_PASSIVE_USER if not set in the Makefile
#ifndef NCOMP_PASSIVE_USER
#  define NCOMP_PASSIVE_USER  0
#endif

// add built-in scalars
#if ( MODEL == HYDRO )

// entropy (or internal energy) for the dual-energy formalism
# ifdef DUAL_ENERGY
#  define NCOMP_PASSIVE_BUILTIN0    1
# else
#  define NCOMP_PASSIVE_BUILTIN0    0
# endif

// cosmic rays
# ifdef COSMIC_RAY
#  define NCOMP_PASSIVE_BUILTIN1    1
# else
#  define NCOMP_PASSIVE_BUILTIN1    0
# endif

// total number of built-in scalars
#  define NCOMP_PASSIVE_BUILTIN     ( NCOMP_PASSIVE_BUILTIN0 + NCOMP_PASSIVE_BUILTIN1 )

#endif // #if ( MODEL == HYDRO )

// define NCOMP_PASSIVE_BUILTIN if not set yet
#ifndef NCOMP_PASSIVE_BUILTIN
#  define NCOMP_PASSIVE_BUILTIN     0
#endif

// total number of passive scalars
#  define NCOMP_PASSIVE       ( NCOMP_PASSIVE_USER + NCOMP_PASSIVE_BUILTIN )

// assuming all passive scalars have the corresponding fluxes
#  define NFLUX_PASSIVE       NCOMP_PASSIVE


// total number of variables in each cell and in the flux array including both active and passive variables
#  define NCOMP_TOTAL         ( NCOMP_FLUID + NCOMP_PASSIVE )
#  define NFLUX_TOTAL         ( NFLUX_FLUID + NFLUX_PASSIVE )


// number of input/output fluid variables in the fluid solver
#if   ( MODEL == HYDRO )
#  define FLU_NIN             NCOMP_TOTAL
#  define FLU_NOUT            NCOMP_TOTAL

// for ELBDM, we do not need to transfer the density component into GPU
#elif ( MODEL == ELBDM )
#  define FLU_NIN             ( NCOMP_TOTAL - 1 )
#  define FLU_NOUT            ( NCOMP_TOTAL - 0 )

#elif ( MODEL == PAR_ONLY )
#  define FLU_NIN             0
#  define FLU_NOUT            0

#else
#  error : ERROR : unsupported MODEL (please edit FLU_NIN and FLU_NOUT for the new MODEL) !!
#endif // MODEL


// number of input fluid variables in the dt solver
// --> EOS_GAMMA/EOS_ISOTHERMAL do not require passive scalars
#if (  MODEL == HYDRO  &&  ( EOS == EOS_GAMMA || EOS == EOS_ISOTHERMAL )  )
#  define FLU_NIN_T           NCOMP_FLUID
#else
#  define FLU_NIN_T           NCOMP_TOTAL
#endif


// number of input/output fluid variables in the source-term solver
// --> fixed to NCOMP_TOTAL for now
#  define FLU_NIN_S           NCOMP_TOTAL
#  define FLU_NOUT_S          NCOMP_TOTAL


// maximum number of output derived fields
#  define DER_NOUT_MAX        10


// built-in fields in different models
#if   ( MODEL == HYDRO )
// field indices of fluid[] --> element of [0 ... NCOMP_FLUID-1]
// --> must NOT modify their values
// --> in addition, they must be consistent with the order these fields are declared in Init_Field()
#  define DENS                0
#  define MOMX                1
#  define MOMY                2
#  define MOMZ                3
#  define ENGY                4

// field indices of passive[] --> element of [NCOMP_FLUID ... NCOMP_TOTAL-1]
#if ( NCOMP_PASSIVE > 0 )

// always put the built-in variables at the END of the field list
// --> so that their indices (e.g., ENPY/EINT/CRAY) can be determined during compilation
// --> convenient (and probably also more efficient) for the fluid solver
#  define PASSIVE_NEXT_IDX0   ( NCOMP_TOTAL - 1   )

# if   ( DUAL_ENERGY == DE_ENPY )
#  define ENPY                ( PASSIVE_NEXT_IDX0 )
#  define PASSIVE_NEXT_IDX1   ( ENPY - 1          )
# elif ( DUAL_ENERGY == DE_EINT )
#  define EINT                ( PASSIVE_NEXT_IDX0 )
#  define PASSIVE_NEXT_IDX1   ( EINT - 1          )
# else
#  define PASSIVE_NEXT_IDX1   ( PASSIVE_NEXT_IDX0 )
# endif

# ifdef COSMIC_RAY
#  define CRAY                ( PASSIVE_NEXT_IDX1 )
#  define PASSIVE_NEXT_IDX2   ( CRAY - 1          )
# else
#  define PASSIVE_NEXT_IDX2   ( PASSIVE_NEXT_IDX1 )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// field indices of magnetic --> element of [0 ... NCOMP_MAG-1]
# ifdef MHD
#  define MAGX                0
#  define MAGY                1
#  define MAGZ                2
# endif

// flux indices of flux[] --> element of [0 ... NFLUX_FLUID-1]
#  define FLUX_DENS           0
#  define FLUX_MOMX           1
#  define FLUX_MOMY           2
#  define FLUX_MOMZ           3
#  define FLUX_ENGY           4

// flux indices of flux_passive[] --> element of [NFLUX_FLUID ... NFLUX_TOTAL-1]
#if ( NCOMP_PASSIVE > 0 )

// always put the built-in variables at the END of the list
#  define FLUX_NEXT_IDX0   ( NFLUX_TOTAL - 1 )

# if   ( DUAL_ENERGY == DE_ENPY )
#  define FLUX_ENPY        ( FLUX_NEXT_IDX0  )
#  define FLUX_NEXT_IDX1   ( FLUX_ENPY - 1   )
# elif ( DUAL_ENERGY == DE_EINT )
#  define FLUX_EINT        ( FLUX_NEXT_IDX0  )
#  define FLUX_NEXT_IDX1   ( FLUX_EINT - 1   )
# else
#  define FLUX_NEXT_IDX1   ( FLUX_NEXT_IDX0  )
# endif

# ifdef COSMIC_RAY
#  define FLUX_CRAY        ( FLUX_NEXT_IDX1  )
#  define FLUX_NEXT_IDX2   ( FLUX_CRAY - 1   )
# else
#  define FLUX_NEXT_IDX2   ( FLUX_NEXT_IDX1  )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// bitwise field indices
// --> must have "_VAR_NAME = 1L<<VAR_NAME" (e.g., _DENS == 1L<<DENS)
// --> convenient for determining subsets of fields (e.g., _DENS|_ENGY)
// --> used as function parameters (e.g., Prepare_PatchData(), Flu_FixUp(), Flu_FixUp_Restrict(), Buf_GetBufferData())
#  define _DENS               ( 1L << DENS )
#  define _MOMX               ( 1L << MOMX )
#  define _MOMY               ( 1L << MOMY )
#  define _MOMZ               ( 1L << MOMZ )
#  define _ENGY               ( 1L << ENGY )

#if ( NCOMP_PASSIVE > 0 )

# if   ( DUAL_ENERGY == DE_ENPY )
#  define _ENPY               ( 1L << ENPY )
# elif ( DUAL_ENERGY == DE_EINT )
#  define _EINT               ( 1L << EINT )
# endif

# ifdef COSMIC_RAY
#  define _CRAY               ( 1L << CRAY )
# endif

#endif // #if ( NCOMP_PASSIVE > 0 )

// magnetic field
# ifdef MHD
#  define _MAGX               ( 1L << MAGX )
#  define _MAGY               ( 1L << MAGY )
#  define _MAGZ               ( 1L << MAGZ )
#  define _MAG                ( _MAGX | _MAGY | _MAGZ )
# else
#  define _MAG                0
# endif

// bitwise flux indices
#  define _FLUX_DENS          ( 1L << FLUX_DENS )
#  define _FLUX_MOMX          ( 1L << FLUX_MOMX )
#  define _FLUX_MOMY          ( 1L << FLUX_MOMY )
#  define _FLUX_MOMZ          ( 1L << FLUX_MOMZ )
#  define _FLUX_ENGY          ( 1L << FLUX_ENGY )

#if ( NFLUX_PASSIVE > 0 )

# if   ( DUAL_ENERGY == DE_ENPY )
#  define _FLUX_ENPY          ( 1L << FLUX_ENPY )
# elif ( DUAL_ENERGY == DE_EINT )
#  define _FLUX_EINT          ( 1L << FLUX_EINT )
# endif

# ifdef COSMIC_RAY
#  define _FLUX_CRAY          ( 1L << FLUX_CRAY )
# endif

#endif // #if ( NFLUX_PASSIVE > 0 )

// bitwise indices of derived fields
// --> start from (1L<<NCOMP_TOTAL) to distinguish from the intrinsic fields
// --> remember to define NDERIVE = total number of derived fields
// _EINT_DER is a derived field for distinguishing from _EINT
// --> the latter is an intrinsic field when adopting DUAL_ENERGY == DE_EINT
#  define _VELX               ( 1L << (NCOMP_TOTAL+ 0) )
#  define _VELY               ( 1L << (NCOMP_TOTAL+ 1) )
#  define _VELZ               ( 1L << (NCOMP_TOTAL+ 2) )
#  define _VELR               ( 1L << (NCOMP_TOTAL+ 3) )
#  define _PRES               ( 1L << (NCOMP_TOTAL+ 4) )
#  define _TEMP               ( 1L << (NCOMP_TOTAL+ 5) )
#  define _EINT_DER           ( 1L << (NCOMP_TOTAL+ 6) )
#  define _MAGX_CC            ( 1L << (NCOMP_TOTAL+ 7) )
#  define _MAGY_CC            ( 1L << (NCOMP_TOTAL+ 8) )
#  define _MAGZ_CC            ( 1L << (NCOMP_TOTAL+ 9) )
#  define _MAG_ENGY_CC        ( 1L << (NCOMP_TOTAL+10) )
#  define _DERIVED            ( _VELX | _VELY | _VELZ | _VELR | _PRES | _TEMP | _EINT_DER | _MAGX_CC | _MAGY_CC | _MAGZ_CC | _MAG_ENGY_CC )
#  define NDERIVE             11


#elif ( MODEL == ELBDM )
// field indices of fluid[] --> element of [0 ... NCOMP_FLUID-1]
#  define  DENS               0
#  define  REAL               1
#  define  IMAG               2

// field indices of passive[] --> element of [NCOMP_FLUID ... NCOMP_TOTAL-1]
// none for ELBDM

// flux indices of flux[] --> element of [0 ... NFLUX_FLUID-1]
#  define  FLUX_DENS          0

// bitwise field indices
#  define _DENS               ( 1L << DENS )
#  define _REAL               ( 1L << REAL )
#  define _IMAG               ( 1L << IMAG )

// bitwise flux indices
#  define _FLUX_DENS          ( 1L << FLUX_DENS )

// bitwise indices of derived fields
#  define _DERIVED            0
#  define NDERIVE             0


#elif ( MODEL == PAR_ONLY )
#  define _DERIVED            0
#  define NDERIVE             0


#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// bitwise field indices used by all models
#  define _NONE               0
# ifdef GRAVITY
#  define _POTE               ( 1L << (NCOMP_TOTAL+NDERIVE) )
# endif
#  define _FLUID              (  ( 1L << NCOMP_FLUID ) - 1L           )
#  define _PASSIVE            (  ( 1L << NCOMP_TOTAL ) - 1L - _FLUID  )
#  define _TOTAL              (  ( 1L << NCOMP_TOTAL ) - 1L           )

#  define _FLUX_FLUID         (  ( 1L << NFLUX_FLUID ) - 1L                )
#  define _FLUX_PASSIVE       (  ( 1L << NFLUX_TOTAL ) - 1L - _FLUX_FLUID  )
#  define _FLUX_TOTAL         (  ( 1L << NFLUX_TOTAL ) - 1L                )



// symbolic constants for particles
#ifdef PARTICLE

// number of built-in particle attributes
// (1) mass, position*3, velocity*3, time, and type
#  define PAR_NATT_BUILTIN0   9

// acceleration*3 when STORE_PAR_ACC is adopted
# ifdef STORE_PAR_ACC
#  define PAR_NATT_BUILTIN1   3
# else
#  define PAR_NATT_BUILTIN1   0
# endif

// particle creation time when STAR_FORMATION is adopted
# ifdef STAR_FORMATION
#  define PAR_NATT_BUILTIN2   1
# else
#  define PAR_NATT_BUILTIN2   0
# endif

// **total** number of built-in particle attributes
#  define PAR_NATT_BUILTIN    ( PAR_NATT_BUILTIN0 + PAR_NATT_BUILTIN1 + PAR_NATT_BUILTIN2 )


// number of particle attributes that we do not want to store on disk (currently time + acceleration*3)
#  define PAR_NATT_UNSTORED   ( 1 + PAR_NATT_BUILTIN1 )
#  define PAR_NATT_STORED     ( PAR_NATT_TOTAL - PAR_NATT_UNSTORED )


// define PAR_NATT_USER if not set in the Makefile
# ifndef PAR_NATT_USER
#  define PAR_NATT_USER       0
# endif


// total number of particle attributes (built-in + user-defined)
#  define PAR_NATT_TOTAL      ( PAR_NATT_BUILTIN + PAR_NATT_USER )


// indices of built-in particle attributes in Par->Attribute[]
// --> must NOT modify their values
#  define  PAR_MASS           0
#  define  PAR_POSX           1
#  define  PAR_POSY           2
#  define  PAR_POSZ           3
#  define  PAR_VELX           4
#  define  PAR_VELY           5
#  define  PAR_VELZ           6
#  define  PAR_TYPE           7

// always put acceleration and time at the END of the particle attribute list
// --> make it easier to discard them when storing data on disk (see Output_DumpData_Total(_HDF5).cpp)
# ifdef STORE_PAR_ACC
#  define  PAR_ACCX           ( PAR_NATT_TOTAL - 4 )
#  define  PAR_ACCY           ( PAR_NATT_TOTAL - 3 )
#  define  PAR_ACCZ           ( PAR_NATT_TOTAL - 2 )
# endif
#  define  PAR_TIME           ( PAR_NATT_TOTAL - 1 )


// bitwise field indices related to particles
// --> note that _POTE = ( 1L << (NCOMP_TOTAL+NDERIVE) )
#  define _PAR_DENS           ( 1L << (NCOMP_TOTAL+NDERIVE+1) )

# if ( MODEL == PAR_ONLY )
#  define _TOTAL_DENS         ( _PAR_DENS )
# else
#  define _TOTAL_DENS         ( 1L << (NCOMP_TOTAL+NDERIVE+2) )
# endif

// particle type macros

// number of particle types MUST be 3 for now
#  define  PAR_NTYPE          3

#  define  PTYPE_TRACER       (real)0
#  define  PTYPE_DARK_MATTER  (real)1
#  define  PTYPE_STAR         (real)2

#else // #ifdef PARTICLE

// total density equals gas density if there is no particle
#  define _TOTAL_DENS         ( _DENS )

#endif // #ifdef PARTICLE ... else ...



// number of fluid ghost zones for the fluid solver
#if   ( MODEL == HYDRO )   // hydro
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
#    if   ( LR_SCHEME == PLM )
#     define LR_GHOST_SIZE          1
#    elif ( LR_SCHEME == PPM )
#     define LR_GHOST_SIZE          2
#    else
#     error : ERROR : unsupported LR_SCHEME !!
#    endif
#  endif // MHM/MHM_RP/CTU

#  if   ( FLU_SCHEME == RTVD )
#     define FLU_GHOST_SIZE         3
#  elif ( FLU_SCHEME == MHM )
#     define FLU_GHOST_SIZE         ( 1 + LR_GHOST_SIZE )
#  elif ( FLU_SCHEME == MHM_RP )
#     define FLU_GHOST_SIZE         ( 2 + LR_GHOST_SIZE )
#  elif ( FLU_SCHEME == CTU )
#    ifdef MHD
#     define FLU_GHOST_SIZE         ( 2 + LR_GHOST_SIZE )
#    else
#     define FLU_GHOST_SIZE         ( 1 + LR_GHOST_SIZE )
#    endif // MHD
#  endif // FLU_SCHEME

#elif ( MODEL == ELBDM )   // ELBDM
#  ifdef LAPLACIAN_4TH
#     define FLU_GHOST_SIZE         6
#  else
#     define FLU_GHOST_SIZE         3
#  endif

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL



// self-gravity constants
#ifdef GRAVITY

// number of input and output variables in the gravity solver
#  if   ( MODEL == HYDRO )
#        define GRA_NIN             NCOMP_FLUID

// for ELBDM, we do not need to transfer the density component
#  elif ( MODEL == ELBDM )
#        define GRA_NIN             ( NCOMP_FLUID - 1 )

#  else
#     error Error : unsupported MODEL (please edit GRA_NIN in the new MODEL) !!
#  endif // MODEL


// number of potential ghost zones for evaluating potential (maximum=5) ~ Poisson solver
#        define POT_GHOST_SIZE      5


// number of potential ghost zones for advancing fluid by gravity ~ Gravity solver
#  if   ( MODEL == HYDRO )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      1
//#      define GRA_GHOST_SIZE      2
#     endif

#  elif ( MODEL == ELBDM )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      0
#     endif

#  elif ( MODEL == PAR_ONLY )
#     ifdef STORE_POT_GHOST
#        define GRA_GHOST_SIZE      2
#     else
#        define GRA_GHOST_SIZE      0
#     endif

#  else
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL


// number of potential ghost zones for correcting the half-step velocity if UNSPLIT_GRAVITY is on
// _F/_G: fluid/gravity solvers
#  ifdef UNSPLIT_GRAVITY
#     if   ( MODEL == HYDRO )
#       ifdef MHD
#        define USG_GHOST_SIZE_F    2
#        define USG_GHOST_SIZE_G    1
#       else
#        define USG_GHOST_SIZE_F    1
#        define USG_GHOST_SIZE_G    1
#       endif
#     elif ( MODEL == ELBDM )
#        define USG_GHOST_SIZE_F    0
#        define USG_GHOST_SIZE_G    0
#     else
#        error : ERROR : unsupported MODEL !!
#     endif // MODEL
#  endif // #ifdef UNSPLIT_GRAVITY


// number of density ghost zones for storing the temporary particle mass density in rho_ext[]
#  ifdef PARTICLE
#        define RHOEXT_GHOST_SIZE   2
#  endif


// number of density ghost zones for the Poisson solver
#        define RHO_GHOST_SIZE      ( POT_GHOST_SIZE-1 )

#endif // #ifdef GRAVITY


// number of ghost zones for the source-term solver
// --> fixed to zero for now since ghost zones in source terms are not supported yet
#        define SRC_GHOST_SIZE      0


// number of ghost zones for computing derived fields
#        define DER_GHOST_SIZE      1


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE                   8
#define PS1             ( 1*PATCH_SIZE )
#define PS2             ( 2*PATCH_SIZE )
#define PS2P1           ( PS2 + 1 )
#define PS1M1           ( PS1 - 1 )
#define PS1P1           ( PS1 + 1 )


// size of GPU arrays (in one dimension)
//###REVISE: support interpolation schemes requiring 2 ghost cells on each side for POT_NXT
#  define FLU_NXT       ( PS2 + 2*FLU_GHOST_SIZE )                // use patch group as the unit
#  define FLU_NXT_P1    ( FLU_NXT + 1 )
#ifdef GRAVITY
#  define POT_NXT       ( PS1/2 + 2*( (POT_GHOST_SIZE+3)/2 ) )    // assuming interpolation ghost zone == 1
#  define RHO_NXT       ( PS1 + 2*RHO_GHOST_SIZE )                // POT/RHO/GRA_NXT use patch as the unit
#  define GRA_NXT       ( PS1 + 2*GRA_GHOST_SIZE )
#  ifdef UNSPLIT_GRAVITY
#  define USG_NXT_F     ( PS2 + 2*USG_GHOST_SIZE_F )              // we use patch group as unit for the fluid   solver
#  define USG_NXT_G     ( PS1 + 2*USG_GHOST_SIZE_G )              // we use patch       as unit for the gravity solver
#  else
#  define USG_NXT_F     ( 1 )                                     // still define USG_NXT_F/G since many function prototypes
#  define USG_NXT_G     ( 1 )                                     // require it
#  endif
#else
#  define GRA_NXT       ( 1 )                                     // still define GRA_NXT   ...
#  define USG_NXT_F     ( 1 )                                     // still define USG_NXT_F ...
#endif
#ifdef PARTICLE
#  define RHOEXT_NXT    ( PS1 + 2*RHOEXT_GHOST_SIZE )             // array rho_ext of each patch
#endif
#  define SRC_NXT       ( PS1 + 2*SRC_GHOST_SIZE )                // use patch as the unit
#  define SRC_NXT_P1    ( SRC_NXT + 1 )
#  define DER_NXT       ( PS1 + 2*DER_GHOST_SIZE )                // use patch as the unit


// size of auxiliary arrays and EoS tables
#if ( MODEL == HYDRO )
#  define EOS_NAUX_MAX           20    // EoS_AuxArray_Flt/Int[]
#  define EOS_NTABLE_MAX         20    // *_EoS_Table[]
#endif

#ifdef GRAVITY
#  define EXT_POT_NAUX_MAX       20    // ExtPot_AuxArray[]
#  define EXT_ACC_NAUX_MAX       20    // ExtAcc_AuxArray[]
#  define EXT_POT_NGENE_MAX       6    // h/d_ExtPotGenePtr
#endif

#if ( MODEL == HYDRO )
#  define SRC_NAUX_DLEP          5     // SrcTerms.Dlep_AuxArray_Flt/Int[]
#  define SRC_DLEP_PROF_NVAR     6     // SrcTerms.Dlep_Profile_DataDevPtr[]/RadiusDevPtr[]
#  define SRC_DLEP_PROF_NBINMAX  4000
#else
#  define SRC_NAUX_DLEP          0
#endif
#  define SRC_NAUX_USER          10    // SrcTerms.User_AuxArray_Flt/Int[]


// bitwise reproducibility in flux and electric field fix-up operations
#if ( MODEL == HYDRO )
# ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_FLUX
# endif

// enable BIT_REP_ELECTRIC by default even when BITWISE_REPRODUCIBILITY is off
// --> ensures that the B field on the common interface between two nearby patches are fully
//     consistent with each other (even the round-off errors are the same)
// --> reduces the div(B) errors significantly
# ifdef MHD
//#ifdef BITWISE_REPRODUCIBILITY
#  define BIT_REP_ELECTRIC
//#endif
# endif // MHD
#endif // HYDRO


// extreme values
#ifndef __INT_MAX__
#  define __INT_MAX__            2147483647
#endif

#ifndef __LONG_MAX__
#  define __LONG_MAX__           9223372036854775807L
#endif

#ifndef __UINT_MAX__
#  define __UINT_MAX__           ( __INT_MAX__*2U + 1U )
#endif

#ifndef __ULONG_MAX__
#  define __ULONG_MAX__          18446744073709551615UL     // 2^64-1
#endif

#ifndef __FLT_MAX__
#  define __FLT_MAX__            3.40282347e+38F
#endif

#ifndef __FLT_MIN__
#  define __FLT_MIN__            1.17549435e-38F
#endif

#ifndef __DBL_MAX__
#  define __DBL_MAX__            1.79769313e+308
#endif

#ifndef __DBL_MIN__
#  define __DBL_MIN__            2.22507386e-308
#endif


// extreme value used for various purposes (e.g., floor value for passive scalars)
#ifdef FLOAT8
#  define TINY_NUMBER            __DBL_MIN__
#  define HUGE_NUMBER            __DBL_MAX__
#else
#  define TINY_NUMBER            __FLT_MIN__
#  define HUGE_NUMBER            __FLT_MAX__
#endif


// maximum allowed error for various purposes (e.g., exact Riemann solver, MHD routines, Mis_CompareRealValue())
#define MAX_ERROR_DBL      1.0e-14
#define MAX_ERROR_FLT      1.0e-06f

#ifdef FLOAT8
#  define MAX_ERROR        MAX_ERROR_DBL
#else
#  define MAX_ERROR        MAX_ERROR_FLT
#endif


// sibling index offset for the non-periodic B.C.
#define SIB_OFFSET_NONPERIODIC   ( -100 )


// son index offset for LOAD_BALANCE
#ifdef LOAD_BALANCE
#  define SON_OFFSET_LB          ( -1000 )
#endif


// flag used in "Buf_RecordBoundaryFlag" and "Flag_Buffer" (must be negative)
#ifndef SERIAL
#  define BUFFER_IS_FLAGGED      ( -999 )
#endif


// marker indicating that the array "pot_ext" has NOT been properly set
#if ( defined GRAVITY  &&  defined STORE_POT_GHOST )
#  define POT_EXT_NEED_INIT      __FLT_MAX__
#endif


// marker indicating that the array "rho_ext" has NOT been properly set
#ifdef PARTICLE
#  define RHO_EXT_NEED_INIT      __FLT_MAX__
#endif


// markers for inactive particles
#ifdef PARTICLE
#  define PAR_INACTIVE_OUTSIDE   ( -1.0 )
#  define PAR_INACTIVE_MPI       ( -2.0 )
#endif


// OpenMP scheduling for particle routines
#if ( defined PARTICLE  &&  defined OPENMP )
#  define PAR_OMP_SCHED          dynamic
#  define PAR_OMP_SCHED_CHUNK    1
#endif


// NULL values
#ifndef NULL
#  define NULL             0
#endif

#ifndef NULL_INT
#  define NULL_INT         __INT_MAX__
#endif

#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__
#endif

#ifndef NULL_BOOL
#  define NULL_BOOL        false
#endif


// GAMER status
#define GAMER_SUCCESS      1
#define GAMER_FAILED       0


// timer switch
#define TIMER_ON           1
#define TIMER_OFF          0


// symbolic constant for Aux_Error()
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// miscellaneous
#define TOP_LEVEL          ( NLEVEL - 1 )


// maximum length for strings
#define MAX_STRING         512




// ############
// ## Macros ##
// ############

// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define   LOG( a )         log( a )
#  define   EXP( a )         exp( a )
#  define  ATAN( a )        atan( a )
#  define FLOOR( a )       floor( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define   POW( a, b )      pow( a, b )
#  define  FMOD( a, b )     fmod( a, b )
#  define ATAN2( a, b )    atan2( a, b )
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define   LOG( a )         logf( a )
#  define   EXP( a )         expf( a )
#  define  ATAN( a )        atanf( a )
#  define FLOOR( a )       floorf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define  FMOD( a, b )     fmodf( a, b )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// square/cube function
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )


// 3D to 1D array indices transformation
#define IDX321( i, j, k, Ni, Nj )   (  ( (k)*(Nj) + (j) )*(Ni) + (i)  )


// 3D to 1D array indices transformation for patch->magnetic[]
#ifdef MHD
#define IDX321_BX( i, j, k, Ni, Nj )   (  ( (k)*((Nj)  ) + (j) )*((Ni)+1) + (i)  )
#define IDX321_BY( i, j, k, Ni, Nj )   (  ( (k)*((Nj)+1) + (j) )*((Ni)  ) + (i)  )
#define IDX321_BZ( i, j, k, Ni, Nj )   (  ( (k)*((Nj)  ) + (j) )*((Ni)  ) + (i)  )
#endif


// helper macros for printing symbolic constants in macros
// ref: https://stackoverflow.com/questions/3419332/c-preprocessor-stringify-the-result-of-a-macro
#  define QUOTE( str )              #str
#  define EXPAND_AND_QUOTE( str )   QUOTE( str )


// convenient macros for defining and declaring global variables
// ==> predefine DEFINE_GLOBAL in the file actually **defines** these global variables
// ==> there should be one and only one file that defines DEFINE_GLOBAL

// SET_GLOBAL will invoke either SET_GLOBAL_INIT or SET_GLOBAL_NOINIT depending on the number of arguments
// ==> http://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments
#define GET_MACRO( _1, _2, TARGET_MACRO, ... )  TARGET_MACRO
#define SET_GLOBAL( ... )   GET_MACRO( __VA_ARGS__, SET_GLOBAL_INIT, SET_GLOBAL_NOINIT ) ( __VA_ARGS__ )

// SET_GLOBAL_INIT/NOINIT are for global variables with/without initialization
#ifdef DEFINE_GLOBAL
# define SET_GLOBAL_INIT( declaration, init_value )   declaration = init_value
# define SET_GLOBAL_NOINIT( declaration )             declaration
#else
# define SET_GLOBAL_INIT( declaration, init_value )   extern declaration
# define SET_GLOBAL_NOINIT( declaration )             extern declaration
#endif


// macro converting an array index (e.g., DENS) to bitwise index (e.g., _DENS=(1L<<DENS))
#define BIDX( idx )     ( 1L << (idx) )



// ################################
// ## Remove useless definitions ##
// ################################
#if ( MODEL == HYDRO )
#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#  undef LR_SCHEME
#  endif

#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#  undef RSOLVER
#  endif
#endif

#if ( MODEL == PAR_ONLY )
#  undef UNSPLIT_GRAVITY
#endif

// currently we always set GPU_ARCH == NONE when GPU is off
#ifndef GPU
#  undef  GPU_ARCH
#  define GPU_ARCH NONE
#endif



#endif  // #ifndef __MACRO_H__
