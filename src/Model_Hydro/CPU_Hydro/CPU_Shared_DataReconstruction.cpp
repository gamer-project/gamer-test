#ifndef __CUFLU_DATARECONSTRUCTION__
#define __CUFLU_DATARECONSTRUCTION__



#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut );
void Hydro_Pri2Con( const real In[], real Out[], const bool FracPassive, const int NFrac, const int FracIdx[],
                    const EoS_DP2E_t EoS_DensPres2Eint, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn );
#if ( FLU_SCHEME == MHM )
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real* const PresIn );
#endif

#endif // #ifdef __CUDACC__ ... else ...


// internal functions (GPU_DEVICE is defined in CUFLU.h)
GPU_DEVICE
static void Hydro_LimitSlope( const real L[], const real C[], const real R[], const LR_Limiter_t LR_Limiter,
                              const real MinMod_Coeff, const int XYZ,
                              const real LEigenVec[][NWAVE], const real REigenVec[][NWAVE], real Slope_Limiter[],
                              const EoS_t *EoS );
#if (  FLU_SCHEME == CTU  ||  ( defined MHD && defined CHAR_RECONSTRUCTION )  )
#ifdef MHD
GPU_DEVICE
static void   MHD_GetEigenSystem( const real CC_Var[], real EigenVal[],
                                  real LEigenVec[][NWAVE], real REigenVec[][NWAVE],
                                  const EoS_t *EoS, const int XYZ );
#else
GPU_DEVICE
static void Hydro_GetEigenSystem( const real CC_Var[], real EigenVal[][NWAVE],
                                  real LEigenVec[][NWAVE], real REigenVec[][NWAVE],
                                  const EoS_t *EoS );
#endif
#endif
#if ( FLU_SCHEME == MHM )
GPU_DEVICE
static void Hydro_HancockPredict( real fc[][NCOMP_LR], const real dt, const real dh,
                                  const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                                  const real MinDens, const real MinPres, const real MinEint,
                                  const EoS_t *EoS );
#endif
#ifdef CHAR_RECONSTRUCTION
GPU_DEVICE
static void Hydro_Pri2Char( real InOut[], const real Dens, const real Pres, const real LEigenVec[][NWAVE],
                            const int XYZ, const EoS_t *EoS );
GPU_DEVICE
static void Hydro_Char2Pri( real InOut[], const real Dens, const real Pres, const real REigenVec[][NWAVE],
                            const int XYZ, const EoS_t *EoS );
#endif


// macro for adding MHD source terms in CTU
#if ( defined MHD  &&  FLU_SCHEME == CTU )
#  define MINMOD( a , b )  (  ( (a)*(b)>(real)0.0 ) ? ( SIGN(a)*FMIN(FABS(a),FABS(b)) ) : (real)0.0  )
#endif




#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. Input data can be either conserved or primitive variables
//                   --> If the input data are conserved variables, one must provide g_ConVar[] and enable "Con2Pri"
//                       --> Primitive variables will be calculated by this function and stored in g_PriVar[]
//                       --> g_PriVar[] must be allocated in advance but it doesn't need to be initialized
//                       --> Adopted by MHM and CTU
//                   --> If the input data are primitive variables, one must provide g_PriVar[] and disable "Con2Pri"
//                       --> g_ConVar[] is useless here
//                       --> Adopted by MHM_RP, where Hydro_RiemannPredict() already returns primitive variables
//                3. Output data are always conserved variables
//                   --> Because Hydro_HancockPredict() only works with conserved variables
//                4. PLM and PPM data reconstruction functions share the same function name
//                5. Face-centered variables will be advanced by half time-step for the MHM and CTU schemes
//                6. Data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION" in the header CUFLU.h
//                7. This function is shared by MHM, MHM_RP, and CTU schemes
//                8. g_FC_B[] has the size of SQR(FLU_NXT)*FLU_NXT_P1 but is accessed with the strides
//                   NIn/NIn+1 along the transverse/longitudinal directions
//                9. Support applying data reconstruction to internal energy and using that instead of pressure
//                   for converting primitive variables to conserved variables
//                   --> Controlled by the option "LR_EINT" in CUFLU.h; see the description thereof for details
//
// Parameter   :  g_ConVar           : Array storing the input cell-centered conserved variables
//                                     --> Should contain NCOMP_TOTAL variables
//                g_FC_B             : Array storing the input face-centered magnetic field (for MHD only)
//                                     --> Should contain NCOMP_MAG variables
//                g_PriVar           : Array storing/to store the cell-centered primitive variables
//                                     --> Should contain NCOMP_LR variables
//                                         --> Store internal energy as the last variable when LR_EINT is on
//                                     --> For MHD, this array currently stores the normal B field as well
//                                     --> For MHM, g_ConVar[] and g_PriVar[] must point to different arrays since
//                                         Hydro_HancockPredict() requires the original g_ConVar[]
//                g_FC_Var           : Array to store the output face-centered conserved variables
//                                     --> Should contain NCOMP_TOTAL_PLUS_MAG variables
//                g_Slope_PPM        : Array to store the x/y/z slopes for the PPM reconstruction
//                                     --> Should contain NCOMP_LR variables
//                                         --> Store internal energy as the last variable when LR_EINT is on
//                                     --> Useless for PLM
//                Con2Pri            : Convert conserved variables in g_ConVar[] to primitive variables and
//                                     store the results in g_PriVar[]
//                NIn                : Size of g_PriVar[] along each direction
//                                     --> Can be smaller than FLU_NXT
//                NGhost             : Number of ghost zones
//                                      --> "NIn-2*NGhost" cells will be computed along each direction
//                                      --> Size of g_FC_Var[] is assumed to be "(NIn-2*NGhost)^3"
//                                      --> The reconstructed data at cell (i,j,k) will be stored in g_FC_Var[]
//                                          with the index "(i-NGhost,j-NGhost,k-NGhost)"
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                dt                 : Time interval to advance solution (for the CTU scheme)
//                dh                 : Cell size
//                MinDens/Pres/Eint  : Density, pressure, and internal energy floors
//                FracPassive        : true --> convert passive scalars to mass fraction during data reconstruction
//                NFrac              : Number of passive scalars for the option "FracPassive"
//                FracIdx            : Target variable indices for the option "FracPassive"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS                : EoS object
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh,
                               const real MinDens, const real MinPres, const real MinEint,
                               const bool FracPassive, const int NFrac, const int FracIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff,
                               const EoS_t *EoS )
{

//### NOTE: temporary solution to the bug in cuda 10.1 and 10.2 that incorrectly overwrites didx_cc[]
#  if   ( FLU_SCHEME == MHM )
   const int NIn    = FLU_NXT;
#  elif ( FLU_SCHEME == MHM_RP )
   const int NIn    = N_HF_VAR;
#  elif ( FLU_SCHEME == CTU )
   const int NIn    = FLU_NXT;
#  else
#  error : ERROR : unsupported FLU_SCHEME !!
#  endif
   const int NGhost = LR_GHOST_SIZE;


// check
#  ifdef GAMER_DEBUG
   if ( NIn - 2*NGhost != N_FC_VAR )
      printf( "ERROR : NIn - 2*NGhost != N_FC_VAR (NIn %d, NGhost %d, N_FC_VAR %d) !!\n",
              NIn, NGhost, N_FC_VAR );

#  if ( defined LR_EINT  &&  FLU_SCHEME == CTU )
#     error : CTU does NOT support LR_EINT !!
#  endif
#  endif // GAMER_DEBUG


   const int  didx_cc[3]     = { 1, NIn, SQR(NIn) };

#  if ( FLU_SCHEME == CTU )
   const real dt_dh2         = (real)0.5*dt/dh;

// index mapping between arrays with size NWAVE and NCOMP_TOTAL_PLUS_MAG/NCOMP_LR
#  ifdef MHD
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4, MAG_OFFSET+1, MAG_OFFSET+2 };
#  else
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };
#  endif
#  endif // #if ( FLU_SCHEME == CTU )

// eigenvalues and eigenvectors
// --> constant components of the left and right eigenvector matrices must be initialized
#  if ( FLU_SCHEME == CTU )
   real EigenVal[3][NWAVE];
#  ifdef MHD
   real REigenVec[NWAVE][NWAVE] = { { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       1.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0 },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 1.0,       0.0,       0.0,       0.0, NULL_REAL,       0.0,       0.0 },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };

#  else
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NWAVE][NWAVE] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #ifdef MHD ... else ...

#  elif ( defined MHD  &&  defined CHAR_RECONSTRUCTION ) // #if ( FLU_SCHEME == CTU )
   real EigenVal[3][NWAVE];
   real REigenVec[NWAVE][NWAVE] = { { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       1.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0 },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 1.0,       0.0,       0.0,       0.0, NULL_REAL,       0.0,       0.0 },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };

#  else // #if ( FLU_SCHEME == CTU ) ... elif ...
   real (*const REigenVec)[NWAVE] = NULL;
   real (*const LEigenVec)[NWAVE] = NULL;
#  endif // #if ( FLU_SCHEME ==  CTU ) ... elif ... else ...


// 0. conserved --> primitive variables
   if ( Con2Pri )
   {
      real ConVar_1Cell[NCOMP_TOTAL_PLUS_MAG], PriVar_1Cell[NCOMP_TOTAL_PLUS_MAG];
#     ifdef LR_EINT
      real Eint;
      real* const EintPtr = &Eint;
#     else
      real* const EintPtr = NULL;
#     endif

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

#        ifdef MHD
//       assuming that g_FC_B[] is accessed with the strides NIn/NIn+1 along the transverse/longitudinal directions
         const int size_ij = SQR( NIn );
         const int i       = idx % NIn;
         const int j       = idx % size_ij / NIn;
         const int k       = idx / size_ij;

         MHD_GetCellCenteredBField( ConVar_1Cell+NCOMP_TOTAL, g_FC_B[0], g_FC_B[1], g_FC_B[2], NIn, NIn, NIn, i, j, k );
#        endif

         Hydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, MinPres, FracPassive, NFrac, FracIdx,
                        JeansMinPres, JeansMinPres_Coeff, EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_PriVar[v][idx] = PriVar_1Cell[v];

#        ifdef LR_EINT
         g_PriVar[NCOMP_TOTAL_PLUS_MAG][idx] = Hydro_CheckMinEint( Eint, MinEint ); // store Eint in the last variable
#        endif
      } // CGPU_LOOP( idx, CUBE(NIn) )

#     ifdef __CUDACC__
      __syncthreads();
#     endif
   } // if ( Con2Pri )


// data reconstruction
   const int N_FC_VAR2 = SQR( N_FC_VAR );
#  ifdef MHD
   const int NIn_p1    = NIn + 1;
   int idx_B[NCOMP_MAG];
#  endif

   CGPU_LOOP( idx_fc, CUBE(N_FC_VAR) )
   {
      const int i_cc   = NGhost + idx_fc%N_FC_VAR;
      const int j_cc   = NGhost + idx_fc%N_FC_VAR2/N_FC_VAR;
      const int k_cc   = NGhost + idx_fc/N_FC_VAR2;
      const int idx_cc = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

#     ifdef MHD
//    assuming that g_FC_B[] is accessed with the strides NIn/NIn+1 along the transverse/longitudinal directions
      idx_B[0] = IDX321( i_cc, j_cc, k_cc, NIn_p1, NIn    );
      idx_B[1] = IDX321( i_cc, j_cc, k_cc, NIn,    NIn_p1 );
      idx_B[2] = IDX321( i_cc, j_cc, k_cc, NIn,    NIn    );
#     endif

//    cc_C/L/R: cell-centered variables of the Central/Left/Right cells
//    fc: face-centered variables of the central cell
      real cc_C[NCOMP_LR], cc_L[NCOMP_LR], cc_R[NCOMP_LR];
      real fc[6][NCOMP_LR], Slope_Limiter[NCOMP_LR];

      for (int v=0; v<NCOMP_LR; v++)   cc_C[v] = g_PriVar[v][idx_cc];


//    1-a. evaluate the eigenvalues and eigenvectors along all three directions for the pure-hydro CTU integrator
#     if ( !defined MHD  &&  FLU_SCHEME == CTU )
      Hydro_GetEigenSystem( cc_C, EigenVal, LEigenVec, REigenVec, EoS );
#     endif


//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
//       1-b. evaluate the eigenvalues and eigenvectors along the target direction for the MHD CTU integrator
#        if (  defined MHD  &&  ( FLU_SCHEME == CTU || defined CHAR_RECONSTRUCTION )  )
         MHD_GetEigenSystem( cc_C, EigenVal[d], LEigenVec, REigenVec, EoS, d );
#        endif


//       2. evaluate the monotonic slope
         const int faceL   = 2*d;      // left and right face indices
         const int faceR   = faceL+1;
         const int idx_ccL = idx_cc - didx_cc[d];
         const int idx_ccR = idx_cc + didx_cc[d];

         for (int v=0; v<NCOMP_LR; v++)
         {
            cc_L[v] = g_PriVar[v][idx_ccL];
            cc_R[v] = g_PriVar[v][idx_ccR];
         }

         Hydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, d,
                           LEigenVec, REigenVec, Slope_Limiter, EoS );


//       3. get the face-centered primitive variables
         for (int v=0; v<NCOMP_LR; v++)
         {
            fc[faceL][v] = cc_C[v] - (real)0.5*Slope_Limiter[v];
            fc[faceR][v] = cc_C[v] + (real)0.5*Slope_Limiter[v];
         }

//       ensure the face-centered variables lie between neighboring cell-centered values
         for (int v=0; v<NCOMP_LR; v++)
         {
            real Min, Max;

            Min = ( cc_C[v] < cc_L[v] ) ? cc_C[v] : cc_L[v];
            Max = ( cc_C[v] > cc_L[v] ) ? cc_C[v] : cc_L[v];
            fc[faceL][v] = ( fc[faceL][v] > Min ) ? fc[faceL][v] : Min;
            fc[faceL][v] = ( fc[faceL][v] < Max ) ? fc[faceL][v] : Max;
            fc[faceR][v] = (real)2.0*cc_C[v] - fc[faceL][v];

            Min = ( cc_C[v] < cc_R[v] ) ? cc_C[v] : cc_R[v];
            Max = ( cc_C[v] > cc_R[v] ) ? cc_C[v] : cc_R[v];
            fc[faceR][v] = ( fc[faceR][v] > Min ) ? fc[faceR][v] : Min;
            fc[faceR][v] = ( fc[faceR][v] < Max ) ? fc[faceR][v] : Max;
            fc[faceL][v] = (real)2.0*cc_C[v] - fc[faceR][v];
         }


//       4. advance the face-centered variables by half time-step for the CTU integrator
#        if ( FLU_SCHEME == CTU )

#        ifdef LR_EINT
#           error : CTU does NOT support LR_EINT !!
#        endif

         real Coeff_L, Coeff_R;
         real Correct_L[NCOMP_LR], Correct_R[NCOMP_LR], dfc[NCOMP_LR];

//       4-1. evaluate the slope (for passive scalars as well)
         for (int v=0; v<NCOMP_LR; v++)   dfc[v] = fc[faceR][v] - fc[faceL][v];


//       4-2. re-order variables for the y/z directions
         Hydro_Rotate3D( dfc, d, true, MAG_OFFSET );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC/HLLD)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )

//       4-2-a1. evaluate the corrections to the left and right face-centered variables

         for (int v=0; v<NWAVE; v++)
         {
            Correct_L[ idx_wave[v] ] = (real)0.0;
            Correct_R[ idx_wave[v] ] = (real)0.0;
         }

#        ifdef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NWAVE; Mode++)
         {
            Coeff_L = (real)0.0;

            for (int v=0; v<NWAVE; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[ idx_wave[v] ];

            Coeff_L *= -dt_dh2*EigenVal[d][Mode];

            for (int v=0; v<NWAVE; v++)   Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
         } // for (int Mode=0; Mode<NWAVE; Mode++)

         for (int v=0; v<NWAVE; v++)   Correct_R[ idx_wave[v] ] = Correct_L[ idx_wave[v] ];

#        else // #ifdef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NWAVE; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NWAVE; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[ idx_wave[v] ];

               Coeff_L *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NWAVE; v++)   Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NWAVE; v++)   Coeff_R += LEigenVec[Mode][v]*dfc[ idx_wave[v] ];

               Coeff_R *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NWAVE; v++)   Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NWAVE; Mode++)

#        endif // #ifdef HLL_INCLUDE_ALL_WAVES ... else ...


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT || ifndef HLL_NO_REF_STATE )

//       4-2-b1. evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][       0 ], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][ NWAVE-1 ], (real)0.0 );

         for (int v=0; v<NWAVE; v++)
         {
            Correct_L[ idx_wave[v] ] = Coeff_L*dfc[ idx_wave[v] ];
            Correct_R[ idx_wave[v] ] = Coeff_R*dfc[ idx_wave[v] ];
         }


//       4-2-b2. evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NWAVE; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NWAVE; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[ idx_wave[v] ];

               Coeff_L *= dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );

               for (int v=0; v<NWAVE; v++)   Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NWAVE; v++)   Coeff_R += LEigenVec[Mode][v]*dfc[ idx_wave[v] ];

               Coeff_R *= dt_dh2*( EigenVal[d][ NWAVE-1 ] - EigenVal[d][Mode] );

               for (int v=0; v<NWAVE; v++)   Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NWAVE; Mode++)

#        endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       4-3. evaluate the corrections to the left and right face-centered passive scalars
//            --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*dfc[v];
            Correct_R[v] = Coeff_R*dfc[v];
         }
#        endif


//       4-4. add the MHD source terms
#        ifdef MHD
         const int t1 = (d+1)%3;    // transverse direction 1
         const int t2 = (d+2)%3;    // transverse direction 2
         real B_nL, B_nR, B_t1L, B_t1R, B_t2L, B_t2R;
         real dB_n, dB_t1, dB_t2, v_t1, v_t2, src_t1, src_t2;

         B_nL   = g_FC_B[d ][ idx_B[d ] ];
         B_t1L  = g_FC_B[t1][ idx_B[t1] ];
         B_t2L  = g_FC_B[t2][ idx_B[t2] ];
         B_nR   = g_FC_B[d ][ idx_B[d ] + didx_cc[d ] ];
         B_t1R  = g_FC_B[t1][ idx_B[t1] + didx_cc[t1] ];
         B_t2R  = g_FC_B[t2][ idx_B[t2] + didx_cc[t2] ];

         dB_n   = B_nR  - B_nL;
         dB_t1  = B_t1R - B_t1L;
         dB_t2  = B_t2R - B_t2L;

         v_t1   = cc_C[ 1 + t1 ];
         v_t2   = cc_C[ 1 + t2 ];

         src_t1 = dt_dh2*v_t1*MINMOD( dB_n, -dB_t1 );
         src_t2 = dt_dh2*v_t2*MINMOD( dB_n, -dB_t2 );

         Correct_L[ MAG_OFFSET + 1 ] += src_t1;
         Correct_R[ MAG_OFFSET + 1 ] += src_t1;
         Correct_L[ MAG_OFFSET + 2 ] += src_t2;
         Correct_R[ MAG_OFFSET + 2 ] += src_t2;
#        endif // #ifdef MHD


//       4-5. evaluate the face-centered variables at the half time-step
         Hydro_Rotate3D( Correct_L, d, false, MAG_OFFSET );
         Hydro_Rotate3D( Correct_R, d, false, MAG_OFFSET );

         for (int v=0; v<NCOMP_LR; v++)
         {
            fc[faceL][v] += Correct_L[v];
            fc[faceR][v] += Correct_R[v];
         }


//       4-6. apply density and pressure floors
         fc[faceL][0] = FMAX( fc[faceL][0], MinDens );
         fc[faceR][0] = FMAX( fc[faceR][0], MinDens );

         fc[faceL][4] = Hydro_CheckMinPres( fc[faceL][4], MinPres );
         fc[faceR][4] = Hydro_CheckMinPres( fc[faceR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         fc[faceL][v] = FMAX( fc[faceL][v], TINY_NUMBER );
         fc[faceR][v] = FMAX( fc[faceR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )


//       5. reset the longitudinal B field to the input face-centered values
//          --> actually no data reconstruction is required for that
//###OPTIMIZARION: do not perform data reconstruction for the longitudinal B field
#        ifdef MHD
#        if ( FLU_SCHEME != CTU )
         const real B_nL = g_FC_B[d][ idx_B[d]              ];
         const real B_nR = g_FC_B[d][ idx_B[d] + didx_cc[d] ];
#        endif
         fc[faceL][ MAG_OFFSET + d ] = B_nL;
         fc[faceR][ MAG_OFFSET + d ] = B_nR;
#        endif // #ifdef MHD


//       6. primitive variables --> conserved variables
//          --> When LR_EINT is on, use the reconstructed internal energy instead of pressure in Hydro_Pri2Con()
//              to skip expensive EoS conversion
         real tmp[NCOMP_LR];  // input and output arrays must not overlap for Pri2Con()
#        ifdef LR_EINT
         real* const EintPtr = tmp + NCOMP_TOTAL_PLUS_MAG;
#        else
         real* const EintPtr = NULL;
#        endif

         for (int v=0; v<NCOMP_LR; v++)   tmp[v] = fc[faceL][v];
         Hydro_Pri2Con( tmp, fc[faceL], FracPassive, NFrac, FracIdx, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

         for (int v=0; v<NCOMP_LR; v++)   tmp[v] = fc[faceR][v];
         Hydro_Pri2Con( tmp, fc[faceR], FracPassive, NFrac, FracIdx, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    7. advance the face-centered variables by half time-step for the MHM integrator
      Hydro_HancockPredict( fc, dt, dh, g_ConVar, idx_cc, MinDens, MinPres, MinEint, EoS );
#     endif


//    8. store the face-centered values to the output array
//       --> use NCOMP_TOTAL_PLUS_MAG instead of LR_EINT since we don't need to store internal energy in g_FC_Var[]
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(N_FC_VAR) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  See the PLM routine
//
// Parameter   :  See the PLM routine
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh,
                               const real MinDens, const real MinPres, const real MinEint,
                               const bool FracPassive, const int NFrac, const int FracIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff,
                               const EoS_t *EoS )
{

//### NOTE: temporary solution to the bug in cuda 10.1 and 10.2 that incorrectly overwrites didx_cc[]
#  if   ( FLU_SCHEME == MHM )
   const int NIn    = FLU_NXT;
#  elif ( FLU_SCHEME == MHM_RP )
   const int NIn    = N_HF_VAR;
#  elif ( FLU_SCHEME == CTU )
   const int NIn    = FLU_NXT;
#  else
#  error : ERROR : unsupported FLU_SCHEME !!
#  endif
   const int NGhost = LR_GHOST_SIZE;


// check
#  ifdef GAMER_DEBUG
   if ( NIn - 2*NGhost != N_FC_VAR )
      printf( "ERROR : NIn - 2*NGhost != N_FC_VAR (NIn %d, NGhost %d, N_FC_VAR %d) !!\n",
              NIn, NGhost, N_FC_VAR );

#  if ( N_SLOPE_PPM != N_FC_VAR + 2 )
#     error : ERROR : N_SLOPE_PPM != N_FC_VAR + 2 !!
#  endif

#  if ( defined LR_EINT  &&  FLU_SCHEME == CTU )
#     error : CTU does NOT support LR_EINT !!
#  endif
#  endif // GAMER_DEBUG


   const int  didx_cc   [3]  = { 1, NIn, SQR(NIn) };
   const int  didx_slope[3]  = { 1, N_SLOPE_PPM, SQR(N_SLOPE_PPM) };

#  if ( FLU_SCHEME == CTU )
   const real dt_dh2         = (real)0.5*dt/dh;

// index mapping between arrays with size NWAVE and NCOMP_TOTAL_PLUS_MAG/NCOMP_LR
#  ifdef MHD
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4, MAG_OFFSET+1, MAG_OFFSET+2 };
#  else
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };
#  endif

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )
#  ifdef HLL_INCLUDE_ALL_WAVES
   const bool HLL_Include_All_Waves = true;
#  else
   const bool HLL_Include_All_Waves = false;
#  endif
#  endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )
#  endif // #if ( FLU_SCHEME == CTU )

// eigenvalues and eigenvectors
// --> constant components of the left and right eigenvector matrices must be initialized
#  if ( FLU_SCHEME == CTU )
   real EigenVal[3][NWAVE];
#  ifdef MHD
   real REigenVec[NWAVE][NWAVE] = { { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       1.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0 },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 1.0,       0.0,       0.0,       0.0, NULL_REAL,       0.0,       0.0 },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };

#  else
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NWAVE][NWAVE] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #ifdef MHD ... else ...

#  elif ( defined MHD  &&  defined CHAR_RECONSTRUCTION ) // #if ( FLU_SCHEME == CTU )
   real EigenVal[3][NWAVE];
   real REigenVec[NWAVE][NWAVE] = { { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       1.0,       0.0,       0.0,       0.0,       0.0,       0.0,       0.0 },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    {       0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 1.0,       0.0,       0.0,       0.0, NULL_REAL,       0.0,       0.0 },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL },
                                    { 0.0,       0.0, NULL_REAL, NULL_REAL,       0.0, NULL_REAL, NULL_REAL },
                                    { 0.0, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL, NULL_REAL } };

#  else // #if ( FLU_SCHEME == CTU ) ... elif ...
   real (*const REigenVec)[NWAVE] = NULL;
   real (*const LEigenVec)[NWAVE] = NULL;
#  endif // #if ( FLU_SCHEME ==  CTU ) ... elif ... else ...


// 0. conserved --> primitive variables
   if ( Con2Pri )
   {
      real ConVar_1Cell[NCOMP_TOTAL_PLUS_MAG], PriVar_1Cell[NCOMP_TOTAL_PLUS_MAG];
#     ifdef LR_EINT
      real Eint;
      real* const EintPtr = &Eint;
#     else
      real* const EintPtr = NULL;
#     endif

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

#        ifdef MHD
//       assuming that g_FC_B[] is accessed with the strides NIn/NIn+1 along the transverse/longitudinal directions
         const int size_ij = SQR( NIn );
         const int i       = idx % NIn;
         const int j       = idx % size_ij / NIn;
         const int k       = idx / size_ij;

         MHD_GetCellCenteredBField( ConVar_1Cell+NCOMP_TOTAL, g_FC_B[0], g_FC_B[1], g_FC_B[2], NIn, NIn, NIn, i, j, k );
#        endif

         Hydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, MinPres, FracPassive, NFrac, FracIdx,
                        JeansMinPres, JeansMinPres_Coeff, EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_PriVar[v][idx] = PriVar_1Cell[v];

#        ifdef LR_EINT
         g_PriVar[NCOMP_TOTAL_PLUS_MAG][idx] = Hydro_CheckMinEint( Eint, MinEint ); // store Eint in the last variable
#        endif
      } // CGPU_LOOP( idx, CUBE(NIn) )

#     ifdef __CUDACC__
      __syncthreads();
#     endif
   } // if ( Con2Pri )


// 1. evaluate the monotonic slope of all cells
   const int N_SLOPE_PPM2 = SQR( N_SLOPE_PPM );
   CGPU_LOOP( idx_slope, CUBE(N_SLOPE_PPM) )
   {
      const int i_cc   = NGhost - 1 + idx_slope%N_SLOPE_PPM;
      const int j_cc   = NGhost - 1 + idx_slope%N_SLOPE_PPM2/N_SLOPE_PPM;
      const int k_cc   = NGhost - 1 + idx_slope/N_SLOPE_PPM2;
      const int idx_cc = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

//    cc_C/L/R: cell-centered variables of the Central/Left/Right cells
      real cc_C[NCOMP_LR], cc_L[NCOMP_LR], cc_R[NCOMP_LR], Slope_Limiter[NCOMP_LR];

      for (int v=0; v<NCOMP_LR; v++)   cc_C[v] = g_PriVar[v][idx_cc];

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         const int idx_ccL = idx_cc - didx_cc[d];
         const int idx_ccR = idx_cc + didx_cc[d];

#        if ( defined MHD  &&  defined CHAR_RECONSTRUCTION )
         MHD_GetEigenSystem( cc_C, EigenVal[d], LEigenVec, REigenVec, EoS, d );
#        endif

         for (int v=0; v<NCOMP_LR; v++)
         {
            cc_L[v] = g_PriVar[v][idx_ccL];
            cc_R[v] = g_PriVar[v][idx_ccR];
         }

         Hydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, d,
                           LEigenVec, REigenVec, Slope_Limiter, EoS );

//       store the results to g_Slope_PPM[]
         for (int v=0; v<NCOMP_LR; v++)   g_Slope_PPM[d][v][idx_slope] = Slope_Limiter[v];

      } // for (int d=0; d<3; d++)
   } // CGPU_LOOP( idx_slope, CUBE(N_SLOPE_PPM) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif


// data reconstruction
   const int N_FC_VAR2 = SQR( N_FC_VAR );
#  ifdef MHD
   const int NIn_p1    = NIn + 1;
   int idx_B[NCOMP_MAG];
#  endif

   CGPU_LOOP( idx_fc, CUBE(N_FC_VAR) )
   {
      const int i_fc      = idx_fc%N_FC_VAR;
      const int j_fc      = idx_fc%N_FC_VAR2/N_FC_VAR;
      const int k_fc      = idx_fc/N_FC_VAR2;

      const int i_cc      = i_fc + NGhost;
      const int j_cc      = j_fc + NGhost;
      const int k_cc      = k_fc + NGhost;
      const int idx_cc    = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

      const int i_slope   = i_fc + 1;   // because N_SLOPE_PPM = N_FC_VAR + 2
      const int j_slope   = j_fc + 1;
      const int k_slope   = k_fc + 1;
      const int idx_slope = IDX321( i_slope, j_slope, k_slope, N_SLOPE_PPM, N_SLOPE_PPM );

#     ifdef MHD
//    assuming that g_FC_B[] is accessed with the strides NIn/NIn+1 along the transverse/longitudinal directions
      idx_B[0] = IDX321( i_cc, j_cc, k_cc, NIn_p1, NIn    );
      idx_B[1] = IDX321( i_cc, j_cc, k_cc, NIn,    NIn_p1 );
      idx_B[2] = IDX321( i_cc, j_cc, k_cc, NIn,    NIn    );
#     endif

 //   cc/fc: cell/face-centered variables; _C_ncomp: central cell with all NCOMP_LR variables
      real cc_C_ncomp[NCOMP_LR], fc[6][NCOMP_LR], dfc[NCOMP_LR], dfc6[NCOMP_LR];

      for (int v=0; v<NCOMP_LR; v++)   cc_C_ncomp[v] = g_PriVar[v][idx_cc];


//    2-a. evaluate the eigenvalues and eigenvectors along all three directions for the pure-hydro CTU integrator
#     if ( !defined MHD  &&  FLU_SCHEME == CTU )
      Hydro_GetEigenSystem( cc_C_ncomp, EigenVal, LEigenVec, REigenVec, EoS );
#     endif


//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
//       2-b. evaluate the eigenvalues and eigenvectors along the target direction for the MHD CTU integrator
#        if (  defined MHD  &&  ( FLU_SCHEME == CTU || defined CHAR_RECONSTRUCTION )  )
         MHD_GetEigenSystem( cc_C_ncomp, EigenVal[d], LEigenVec, REigenVec, EoS, d );
#        endif


//       3. get the face-centered primitive variables
         const int faceL      = 2*d;      // left and right face indices
         const int faceR      = faceL+1;
         const int idx_ccL    = idx_cc - didx_cc[d];
         const int idx_ccR    = idx_cc + didx_cc[d];
         const int idx_slopeL = idx_slope - didx_slope[d];
         const int idx_slopeR = idx_slope + didx_slope[d];

         for (int v=0; v<NCOMP_LR; v++)
         {
//          cc/fc: cell/face-centered variables; _C/L/R: Central/Left/Right cells
            real cc_C, cc_L, cc_R, dcc_L, dcc_R, dcc_C, fc_L, fc_R, Max, Min;

//          3-1. parabolic interpolation
            cc_L  = g_PriVar[v][idx_ccL];
            cc_R  = g_PriVar[v][idx_ccR];
            cc_C  = cc_C_ncomp[v];

            dcc_L = g_Slope_PPM[d][v][idx_slopeL];
            dcc_R = g_Slope_PPM[d][v][idx_slopeR];
            dcc_C = g_Slope_PPM[d][v][idx_slope ];

            fc_L  = (real)0.5*( cc_C + cc_L ) - (real)1.0/(real)6.0*( dcc_C - dcc_L );
            fc_R  = (real)0.5*( cc_C + cc_R ) - (real)1.0/(real)6.0*( dcc_R - dcc_C );


//          3-2. monotonicity constraint
//          extra monotonicity check for the CENTRAL limiter since it's not TVD
            if ( LR_Limiter == LR_LIMITER_CENTRAL )
            {
               if ( (cc_C-fc_L)*(fc_L-cc_L) < (real)0.0 )   fc_L = (real)0.5*( cc_C + cc_L );
               if ( (cc_R-fc_R)*(fc_R-cc_C) < (real)0.0 )   fc_R = (real)0.5*( cc_C + cc_R );
            }

            dfc [v] = fc_R - fc_L;
            dfc6[v] = (real)6.0*(  cc_C - (real)0.5*( fc_L + fc_R )  );

            if (  ( fc_R - cc_C )*( cc_C - fc_L ) <= (real)0.0  )
            {
               fc_L = cc_C;
               fc_R = cc_C;
            }
            else if ( dfc[v]*dfc6[v] > +dfc[v]*dfc[v] )
               fc_L = (real)3.0*cc_C - (real)2.0*fc_R;
            else if ( dfc[v]*dfc6[v] < -dfc[v]*dfc[v] )
               fc_R = (real)3.0*cc_C - (real)2.0*fc_L;


//          3-3. ensure the face-centered variables lie between neighboring cell-centered values
            Min  = ( cc_C < cc_L ) ? cc_C : cc_L;
            Max  = ( cc_C > cc_L ) ? cc_C : cc_L;
            fc_L = ( fc_L > Min  ) ? fc_L : Min;
            fc_L = ( fc_L < Max  ) ? fc_L : Max;

            Min  = ( cc_C < cc_R ) ? cc_C : cc_R;
            Max  = ( cc_C > cc_R ) ? cc_C : cc_R;
            fc_R = ( fc_R > Min  ) ? fc_R : Min;
            fc_R = ( fc_R < Max  ) ? fc_R : Max;

            fc[faceL][v] = fc_L;
            fc[faceR][v] = fc_R;

         } // for (int v=0; v<NCOMP_LR; v++)


//       4. advance the face-centered variables by half time-step for the CTU integrator
#        if ( FLU_SCHEME == CTU )

#        ifdef LR_EINT
#           error : CTU does NOT support LR_EINT !!
#        endif

         real Coeff_L, Coeff_R;
         real Correct_L[NCOMP_LR], Correct_R[NCOMP_LR];

//       4-1. compute the PPM coefficient (for the passive scalars as well)
         for (int v=0; v<NCOMP_LR; v++)
         {
            dfc [v] = fc[faceR][v] - fc[faceL][v];
            dfc6[v] = (real)6.0*(  cc_C_ncomp[v] - (real)0.5*( fc[faceL][v] + fc[faceR][v] )  );
         }

//       4-2. re-order variables for the y/z directions
         Hydro_Rotate3D( dfc,  d, true, MAG_OFFSET );
         Hydro_Rotate3D( dfc6, d, true, MAG_OFFSET );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC/HLLD)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )

//       4-2-a1. evaluate the corrections to the left and right face-centered variables
         for (int v=0; v<NWAVE; v++)
         {
            Correct_L[ idx_wave[v] ] = (real)0.0;
            Correct_R[ idx_wave[v] ] = (real)0.0;
         }

         for (int Mode=0; Mode<NWAVE; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] <= (real)0.0 )
            {
               const real Coeff_C = -dt_dh2*EigenVal[d][Mode];
               const real Coeff_D = real(-4.0/3.0)*SQR(Coeff_C);

               for (int v=0; v<NWAVE; v++)
                  Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[ idx_wave[v] ] + dfc6[ idx_wave[v] ] ) +
                                                   Coeff_D*( dfc6[ idx_wave[v] ]                      )  );

               for (int v=0; v<NWAVE; v++)
                  Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
            }

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] >= (real)0.0 )
            {
               const real Coeff_A = -dt_dh2*EigenVal[d][Mode];
               const real Coeff_B = real(-4.0/3.0)*SQR(Coeff_A);

               for (int v=0; v<NWAVE; v++)
                  Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[ idx_wave[v] ] - dfc6[ idx_wave[v] ] ) +
                                                   Coeff_B*( dfc6[ idx_wave[v] ]                      )  );

               for (int v=0; v<NWAVE; v++)
                  Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NWAVE; Mode++)


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT || ifndef HLL_NO_REF_STATE )

//       4-2-b1. evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][       0 ], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][ NWAVE-1 ], (real)0.0 );

         for (int v=0; v<NWAVE; v++)
         {
            Correct_L[ idx_wave[v] ] = Coeff_L*(  dfc[ idx_wave[v] ] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[ idx_wave[v] ]  );
            Correct_R[ idx_wave[v] ] = Coeff_R*(  dfc[ idx_wave[v] ] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[ idx_wave[v] ]  );
         }


//       4-2-b2. evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NWAVE; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               const real Coeff_C = dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_D=0 when Coeff_C=0
//             Coeff_D = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][   0]*EigenVal[d][   0] -
//                                                      EigenVal[d][Mode]*EigenVal[d][Mode]   );
               const real Coeff_D = real(4.0/3.0)*dt_dh2*Coeff_C*( EigenVal[d][0] + EigenVal[d][Mode] );

               for (int v=0; v<NWAVE; v++)
                  Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[ idx_wave[v] ] + dfc6[ idx_wave[v] ] ) +
                                                   Coeff_D*( dfc6[ idx_wave[v] ]                      )  );

               for (int v=0; v<NWAVE; v++)
                  Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               const real Coeff_A = dt_dh2*( EigenVal[d][ NWAVE-1 ] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_B=0 when Coeff_A=0
//             Coeff_B = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][NWAVE-1]*EigenVal[d][NWAVE-1] -
//                                                      EigenVal[d][Mode   ]*EigenVal[d][Mode   ]   );
               const real Coeff_B = real(4.0/3.0)*dt_dh2*Coeff_A*( EigenVal[d][ NWAVE-1 ] + EigenVal[d][Mode] );

               for (int v=0; v<NWAVE; v++)
                  Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[ idx_wave[v] ] - dfc6[ idx_wave[v] ] ) +
                                                   Coeff_B*( dfc6[ idx_wave[v] ]                      )  );

               for (int v=0; v<NWAVE; v++)
                  Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NWAVE; Mode++)

#        endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       4-3. evaluate the corrections to the left and right face-centered passive scalars
//            --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*(  dfc[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[v]  );
            Correct_R[v] = Coeff_R*(  dfc[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[v]  );
         }
#        endif


//       4-4. add the MHD source terms
#        ifdef MHD
         const int t1 = (d+1)%3;    // transverse direction 1
         const int t2 = (d+2)%3;    // transverse direction 2
         real B_nL, B_nR, B_t1L, B_t1R, B_t2L, B_t2R;
         real dB_n, dB_t1, dB_t2, v_t1, v_t2, src_t1, src_t2;

         B_nL   = g_FC_B[d ][ idx_B[d ] ];
         B_t1L  = g_FC_B[t1][ idx_B[t1] ];
         B_t2L  = g_FC_B[t2][ idx_B[t2] ];
         B_nR   = g_FC_B[d ][ idx_B[d ] + didx_cc[d ] ];
         B_t1R  = g_FC_B[t1][ idx_B[t1] + didx_cc[t1] ];
         B_t2R  = g_FC_B[t2][ idx_B[t2] + didx_cc[t2] ];

         dB_n   = B_nR  - B_nL;
         dB_t1  = B_t1R - B_t1L;
         dB_t2  = B_t2R - B_t2L;

         v_t1   = cc_C_ncomp[ 1 + t1 ];
         v_t2   = cc_C_ncomp[ 1 + t2 ];

         src_t1 = dt_dh2*v_t1*MINMOD( dB_n, -dB_t1 );
         src_t2 = dt_dh2*v_t2*MINMOD( dB_n, -dB_t2 );

         Correct_L[ MAG_OFFSET + 1 ] += src_t1;
         Correct_R[ MAG_OFFSET + 1 ] += src_t1;
         Correct_L[ MAG_OFFSET + 2 ] += src_t2;
         Correct_R[ MAG_OFFSET + 2 ] += src_t2;
#        endif // #ifdef MHD


//       4-5. evaluate the face-centered variables at the half time-step
         Hydro_Rotate3D( Correct_L, d, false, MAG_OFFSET );
         Hydro_Rotate3D( Correct_R, d, false, MAG_OFFSET );

         for (int v=0; v<NCOMP_LR; v++)
         {
            fc[faceL][v] += Correct_L[v];
            fc[faceR][v] += Correct_R[v];
         }


//       4-6. apply density and pressure floors
         fc[faceL][0] = FMAX( fc[faceL][0], MinDens );
         fc[faceR][0] = FMAX( fc[faceR][0], MinDens );

         fc[faceL][4] = Hydro_CheckMinPres( fc[faceL][4], MinPres );
         fc[faceR][4] = Hydro_CheckMinPres( fc[faceR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         fc[faceL][v] = FMAX( fc[faceL][v], TINY_NUMBER );
         fc[faceR][v] = FMAX( fc[faceR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )


//       5. reset the longitudinal B field to the input face-centered values
//          --> actually no data reconstruction is required for that
//###OPTIMIZARION: do not perform data reconstruction for the longitudinal B field
#        ifdef MHD
#        if ( FLU_SCHEME != CTU )
         const real B_nL = g_FC_B[d][ idx_B[d]              ];
         const real B_nR = g_FC_B[d][ idx_B[d] + didx_cc[d] ];
#        endif
         fc[faceL][ MAG_OFFSET + d ] = B_nL;
         fc[faceR][ MAG_OFFSET + d ] = B_nR;
#        endif // #ifdef MHD


//       6. primitive variables --> conserved variables
//          --> When LR_EINT is on, use the reconstructed internal energy instead of pressure in Hydro_Pri2Con()
//              to skip expensive EoS conversion
         real tmp[NCOMP_LR];  // input and output arrays must not overlap for Pri2Con()
#        ifdef LR_EINT
         real* const EintPtr = tmp + NCOMP_TOTAL_PLUS_MAG;
#        else
         real* const EintPtr = NULL;
#        endif

         for (int v=0; v<NCOMP_LR; v++)   tmp[v] = fc[faceL][v];
         Hydro_Pri2Con( tmp, fc[faceL], FracPassive, NFrac, FracIdx, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

         for (int v=0; v<NCOMP_LR; v++)   tmp[v] = fc[faceR][v];
         Hydro_Pri2Con( tmp, fc[faceR], FracPassive, NFrac, FracIdx, EoS->DensPres2Eint_FuncPtr,
                        EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, EintPtr );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    7. advance the face-centered variables by half time-step for the MHM integrator
      Hydro_HancockPredict( fc, dt, dh, g_ConVar, idx_cc, MinDens, MinPres, MinEint, EoS );
#     endif


//    8. store the face-centered values to the output array
//       --> use NCOMP_TOTAL_PLUS_MAG instead of LR_EINT since we don't need to store internal energy in g_FC_Var[]
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(N_FC_VAR) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )



#ifdef CHAR_RECONSTRUCTION
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Pri2Char
// Description :  Primitive variables --> characteristic variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//                2. Input and output share the same array
//                3. InOut[] should have the size of NCOMP_TOTAL_PLUS_MAG or NCOMP_EINT
//                   --> For LR_EINT, where NCOMP_EINT=NCOMP_TOTAL_PLUS_MAG+1, this function assumes that the
//                       internal energy is stored as the last element and does not touch it at all
//                4. Does NOT support general EoS
//
// Parameter   :  InOut     : Array storing both the input primitive variables and output characteristic variables
//                Dens      : Density
//                Pres      : Pressure
//                LEigenVec : Left eigenvector (for MHD only)
//                XYZ       : Target spatial direction : (0/1/2) --> (x/y/z)
//                EoS       : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Pri2Char( real InOut[], const real Dens, const real Pres, const real LEigenVec[][NWAVE],
                     const int XYZ, const EoS_t *EoS )
{

// check
#  if ( EOS == EOS_GAMMA )
#  ifndef MHD
   const real *Passive = NULL;   // EOS_GAMMA does not involve passive scalars
#  endif
#  else
#  error : Hydro_Pri2Char() only supports EOS_GAMMA !!
#  endif

#  if ( defined CHECK_UNPHYSICAL_IN_FLUID  &&  !defined MHD )
   Hydro_CheckUnphysical( NULL, NULL, &Pres, NULL, "pressure", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
#  endif


// back-up the input array and rotate it according to the target direction
// --> it's unnecessary to copy the passive scalars since they will not be modified
   real Temp[ NCOMP_FLUID + NCOMP_MAG ];

   for (int v=0; v<NCOMP_FLUID; v++)   Temp[ v ]               = InOut[ v ];
#  ifdef MHD
   for (int v=0; v<NCOMP_MAG;   v++)   Temp[ v + NCOMP_FLUID ] = InOut[ v + MAG_OFFSET ];
#  endif

   Hydro_Rotate3D( Temp, XYZ, true, NCOMP_FLUID );

// remove the normal B field to be consistent with the eigenvector matrix
#  ifdef MHD
   Temp[ NCOMP_FLUID + 0 ] = Temp[ NCOMP_FLUID + 1 ];
   Temp[ NCOMP_FLUID + 1 ] = Temp[ NCOMP_FLUID + 2 ];
#  endif


// primitive --> characteristic
// a. MHD
#  ifdef MHD
   const real tmp_f1 = LEigenVec[0][1]*Temp[1] + LEigenVec[0][2]*Temp[2] + LEigenVec[0][3]*Temp[3];
   const real tmp_b1 = LEigenVec[0][4]*Temp[4] + LEigenVec[0][5]*Temp[5] + LEigenVec[0][6]*Temp[6];
   const real tmp_f2 = LEigenVec[2][1]*Temp[1] + LEigenVec[2][2]*Temp[2] + LEigenVec[2][3]*Temp[3];
   const real tmp_b2 = LEigenVec[2][4]*Temp[4] + LEigenVec[2][5]*Temp[5] + LEigenVec[2][6]*Temp[6];

   InOut[MAG_OFFSET+0] = (real)0.0;
   InOut[           3] = Temp[0] + LEigenVec[3][4]*Temp[4];
   InOut[           1] = LEigenVec[1][2]*Temp[2] + LEigenVec[1][3]*Temp[3] + LEigenVec[1][5]*Temp[5] + LEigenVec[1][6]*Temp[6];
   InOut[MAG_OFFSET+1] = LEigenVec[5][2]*Temp[2] + LEigenVec[5][3]*Temp[3] + LEigenVec[5][5]*Temp[5] + LEigenVec[5][6]*Temp[6];
   InOut[           0] =  tmp_f1 + tmp_b1;
   InOut[           2] =  tmp_f2 + tmp_b2;
   InOut[           4] = -tmp_f2 + tmp_b2;
   InOut[MAG_OFFSET+2] = -tmp_f1 + tmp_b1;

// b. pure hydro
#  else // #ifdef MHD
   const real  a2 = EoS->DensPres2CSqr_FuncPtr( Dens, Pres, Passive, EoS->AuxArrayDevPtr_Flt,
                                                EoS->AuxArrayDevPtr_Int, EoS->Table, NULL );
   const real _a2 = (real)1.0 / a2;
   const real _a  = SQRT( _a2 );

   InOut[0] = -(real)0.5*Dens*_a*Temp[1] + (real)0.5*_a2*Temp[4];
   InOut[1] = Temp[0] - _a2*Temp[4];
   InOut[2] = Temp[2];
   InOut[3] = Temp[3];
   InOut[4] = +(real)0.5*Dens*_a*Temp[1] + (real)0.5*_a2*Temp[4];
#  endif // #ifdef MHD ... else ...

} // FUNCTION : Hydro_Pri2Char



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Char2Pri
// Description :  Characteristic variables --> primitive variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//                2. Input and output share the same array
//                3. InOut[] should have the size of NCOMP_TOTAL_PLUS_MAG or NCOMP_EINT
//                   --> For LR_EINT, where NCOMP_EINT=NCOMP_TOTAL_PLUS_MAG+1, this function assumes that the
//                       internal energy is stored as the last element and does not touch it at all
//                4. Does NOT support general EoS
//
// Parameter   :  InOut     : Array storing both the input characteristic variables and output primitive variables
//                Dens      : Density
//                Pres      : Pressure
//                REigenVec : Right eigenvector (for MHD only)
//                XYZ       : Target spatial direction : (0/1/2) --> (x/y/z)
//                EoS       : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Char2Pri( real InOut[], const real Dens, const real Pres, const real REigenVec[][NWAVE],
                     const int XYZ, const EoS_t *EoS )
{

// check
#  if ( EOS == EOS_GAMMA )
   const real *Passive = NULL;   // EOS_GAMMA does not involve passive scalars
#  else
#  error : Hydro_Char2Pri() only supports EOS_GAMMA !!
#  endif

#  if ( defined CHECK_UNPHYSICAL_IN_FLUID  &&  !defined MHD )
   Hydro_CheckUnphysical( NULL, NULL, &Pres, NULL, "pressure", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &Dens, NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
#  endif


// back-up the input array and rotate it according to the target direction
// --> it's unnecessary to copy the passive scalars since they will not be modified
// --> it's also unnecessary to copy the normal B field (just to be consistent with the eigenvector matrix)
   real Temp[NWAVE];

   for (int v=0; v<NCOMP_FLUID; v++)      Temp[v] = InOut[v];

#  ifdef MHD
   for (int v=NCOMP_FLUID; v<NWAVE; v++)  Temp[v] = InOut[ v - NCOMP_FLUID + MAG_OFFSET + 1 ];
#  endif


// primitive --> characteristic
   const real a2 = EoS->DensPres2CSqr_FuncPtr( Dens, Pres, Passive, EoS->AuxArrayDevPtr_Flt,
                                               EoS->AuxArrayDevPtr_Int, EoS->Table, NULL );

// a. MHD
#  ifdef MHD
   InOut[           0] = REigenVec[0][0]*Temp[0] + REigenVec[2][0]*Temp[2] + Temp[3] +
                         REigenVec[4][0]*Temp[4] + REigenVec[6][0]*Temp[6];
   InOut[           1] = REigenVec[0][1]*Temp[0] + REigenVec[2][1]*Temp[2] + REigenVec[4][1]*Temp[4] +
                         REigenVec[6][1]*Temp[6];
   InOut[           2] = REigenVec[0][2]*Temp[0] + REigenVec[1][2]*Temp[1] + REigenVec[2][2]*Temp[2] +
                         REigenVec[4][2]*Temp[4] + REigenVec[5][2]*Temp[5] + REigenVec[6][2]*Temp[6];
   InOut[           3] = REigenVec[0][3]*Temp[0] + REigenVec[1][3]*Temp[1] + REigenVec[2][3]*Temp[2] +
                         REigenVec[4][3]*Temp[4] + REigenVec[5][3]*Temp[5] + REigenVec[6][3]*Temp[6];
   InOut[           4] = ( InOut[0] - Temp[3] )*a2;
   InOut[MAG_OFFSET+0] = (real)0.0;
   InOut[MAG_OFFSET+1] = REigenVec[0][5]*Temp[0] + REigenVec[1][5]*Temp[1] + REigenVec[2][5]*Temp[2] +
                         REigenVec[4][5]*Temp[4] + REigenVec[5][5]*Temp[5] + REigenVec[6][5]*Temp[6];
   InOut[MAG_OFFSET+2] = REigenVec[0][6]*Temp[0] + REigenVec[1][6]*Temp[1] + REigenVec[2][6]*Temp[2] +
                         REigenVec[4][6]*Temp[4] + REigenVec[5][6]*Temp[5] + REigenVec[6][6]*Temp[6];

// b. pure hydro
#  else // #ifdef MHD
   const real a = SQRT( a2 );

   InOut[0] = Temp[0] + Temp[1] + Temp[4];
   InOut[1] = a/Dens*( -Temp[0] + Temp[4] );
   InOut[2] = Temp[2];
   InOut[3] = Temp[3];
   InOut[4] = a2*( Temp[0] + Temp[4] );
#  endif // #ifdef MHD ... else ...

   Hydro_Rotate3D( InOut, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_Char2Pri
#endif



#if (  FLU_SCHEME == CTU  ||  ( defined MHD && defined CHAR_RECONSTRUCTION )  )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetEigenSystem
// Description :  Evaluate the eigenvalues and left/right eigenvectors
//
// Note        :  1. Input data must be primitive variables
//                2. Constant components of eigenvectors must be set in advance
//                3. Work for the CTU scheme and the characteristic data reconstruction in MHD
//                4. Do not need to consider passive scalars
//                   --> Their eigenmatrices are just identity matrix
//                5. For pure hydro, this function computes the eigenvalues and eigenvectors
//                   along all three spatial directions at once
//                   --> Because eigenvectors along different directions are the same for pure hydro
//                   But for MHD, this function only computes the eigenvalues and eigenvectors
//                   along the spatial direction specified by XYZ
//                   --> Because eigenvectors along different directions are different for MHD
//                6. Does NOT support general EoS
//
// Parameter   :  CC_Var      : Array storing the input cell-centered primitive variables
//                EigenVal    : Array to store the output eigenvalues
//                              --> Hydro: along all three spatial directions
//                                  MHD  : only along the target spatial direction
//                L/REigenVec : Array to store the output left/right eigenvectors
//                EoS         : EoS object
//                XYZ         : Target spatial direction (for MHD only)
//
// Return      :  EigenVal[], L/REigenVec[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
#ifdef MHD
void   MHD_GetEigenSystem( const real CC_Var[], real EigenVal[],
                           real LEigenVec[][NWAVE], real REigenVec[][NWAVE],
                           const EoS_t *EoS, const int XYZ )
#else
void Hydro_GetEigenSystem( const real CC_Var[], real EigenVal[][NWAVE],
                           real LEigenVec[][NWAVE], real REigenVec[][NWAVE],
                           const EoS_t *EoS )
#endif
{

#  if ( EOS == EOS_GAMMA )
   const real *Passive = NULL;   // EOS_GAMMA does not involve passive scalars
#  else
#  error : Hydro/MHD_GetEigenSystem() only supports EOS_GAMMA !!
#  endif

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_CheckUnphysical( NULL, NULL, &CC_Var[4], NULL, "pressure", __FILE__, __FUNCTION__, __LINE__, true );
   Hydro_CheckUnphysical( NULL, NULL, &CC_Var[0], NULL, "density", __FILE__, __FUNCTION__, __LINE__, true );
#  endif


   const real  Rho = CC_Var[0];
   const real _Rho = (real)1.0/Rho;
   const real  a2  = EoS->DensPres2CSqr_FuncPtr( Rho, CC_Var[4], Passive, EoS->AuxArrayDevPtr_Flt,
                                                 EoS->AuxArrayDevPtr_Int, EoS->Table, NULL );
   const real  a   = SQRT( a2 );
   const real _a   = (real)1.0/a;
   const real _a2  = _a*_a;

// a. MHD
#  ifdef MHD
   real Cf2, Cs2, Cf, Cs;
   real PriVar[NCOMP_TOTAL_PLUS_MAG];

   for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)  PriVar[v] = CC_Var[v];

   Hydro_Rotate3D( PriVar, XYZ, true, MAG_OFFSET );

   const real Bx          = PriVar[ MAG_OFFSET + 0 ];
   const real By          = PriVar[ MAG_OFFSET + 1 ];
   const real Bz          = PriVar[ MAG_OFFSET + 2 ];
   const real Bn2         = SQR( By ) + SQR( Bz );
   const real Bn          = SQRT( Bn2 );
   const real Cax2        = SQR( Bx )*_Rho;
   const real Cax         = SQRT( Cax2 );
   const real Cat2        = Bn2*_Rho;
   const real tsum        = Cax2 + Cat2 + a2;                        // Ca^2 + a^2
   const real tdif        = Cax2 + Cat2 - a2;                        // Ca^2 - a^2
   const real Cf2_min_Cs2 = SQRT( SQR(tdif) + (real)4.0*a2*Cat2 );   // Cf^2 - Cs^2

// evaluate the fast/slow wave speed (Cf/Cs)
   if ( Cat2 == (real)0.0 )
   {
      if ( Cax2 == a2 ) {
         Cf2 = a2;
         Cs2 = a2;
      }
      else if ( Cax2 > a2 ) {
         Cf2 = Cax2;
         Cs2 = a2;
      }
      else {
         Cf2 = a2;
         Cs2 = Cax2;
      }
   }

   else
   {
      if ( Cax2 == (real)0.0 ) {
         Cf2 = a2 + Cat2;
         Cs2 = (real)0.0;
      }
      else {
         Cf2 = (real)0.5*( tsum + Cf2_min_Cs2 );
         Cs2 = a2*Cax2/Cf2;   // do not use "Cf2 - Cf2_min_Cs2" to avoid negative values caused by round-off errors
//       Cs2 = Cf2 - Cf2_min_Cs2;
      }
   } // if ( Cat2 == (real)0.0 ) ... else ...

   Cf = SQRT( Cf2 );
   Cs = SQRT( Cs2 );


// eigenvalues along the target spatial direction
   EigenVal[0] = PriVar[1] - Cf;
   EigenVal[1] = PriVar[1] - Cax;
   EigenVal[2] = PriVar[1] - Cs;
   EigenVal[3] = PriVar[1];
   EigenVal[4] = PriVar[1] + Cs;
   EigenVal[5] = PriVar[1] + Cax;
   EigenVal[6] = PriVar[1] + Cf;


// right eigenvectors (rows instead of columns of the matrix REigenVec for better performance)
   const real S          = SIGN( Bx );
   const real sqrt_Rho   = SQRT( Rho );
   const real _sqrt_Rho  = (real)1.0 / sqrt_Rho;
   const real a2_min_Cs2 = a2 - Cs2;
   const real Cf2_min_a2 = Cf2 - a2;

   real beta_y, beta_z, alpha_f, alpha_s;

   if ( Bn == (real)0.0 ) {
      beta_y = (real)1.0;
      beta_z = (real)0.0;
   }
   else {
      const real _Bn = (real)1.0 / Bn;
      beta_y = By * _Bn;
      beta_z = Bz * _Bn;
   }

   if ( Cf2_min_Cs2 == (real)0.0 ) {
      alpha_f = (real)1.0;
      alpha_s = (real)0.0;
   }
   else if ( a2_min_Cs2 <= (real)0.0 ) {
      alpha_f = (real)0.0;
      alpha_s = (real)1.0;
   }
   else if ( Cf2_min_a2 <= (real)0.0 ) {
      alpha_f = (real)1.0;
      alpha_s = (real)0.0;
   }
   else {
#     ifdef CHECK_UNPHYSICAL_IN_FLUID
      Hydro_CheckUnphysical( NULL, NULL, &a2_min_Cs2, NULL, "a2_min_Cs2", __FILE__, __FUNCTION__, __LINE__, true );
      Hydro_CheckUnphysical( NULL, NULL, &Cf2_min_a2, NULL, "Cf2_min_a2", __FILE__, __FUNCTION__, __LINE__, true );
#     endif

      const real _Cf2_min_Cs2 = (real)1.0 / Cf2_min_Cs2;
      alpha_f = SQRT( a2_min_Cs2*_Cf2_min_Cs2 );
      alpha_s = SQRT( Cf2_min_a2*_Cf2_min_Cs2 );
   }

   const real Af         = a * alpha_f * sqrt_Rho;
   const real As         = a * alpha_s * sqrt_Rho;
   const real Cff        = Cf * alpha_f;
   const real Css        = Cs * alpha_s;
   const real Qf         = Cff * S;
   const real Qs         = Css * S;
   const real S_sqrt_Rho = S * sqrt_Rho;

   REigenVec[0][0] =  Rho * alpha_f;
   REigenVec[0][1] = -Cff;
   REigenVec[0][2] =  Qs * beta_y;
   REigenVec[0][3] =  Qs * beta_z;
   REigenVec[0][4] =  REigenVec[0][0] * a2;
   REigenVec[0][5] =  As * beta_y;
   REigenVec[0][6] =  As * beta_z;

   REigenVec[1][2] = -beta_z;
   REigenVec[1][3] =  beta_y;
   REigenVec[1][5] = -S_sqrt_Rho * beta_z;
   REigenVec[1][6] =  S_sqrt_Rho * beta_y;

   REigenVec[2][0] =  Rho * alpha_s;
   REigenVec[2][1] = -Css;
   REigenVec[2][2] = -Qf * beta_y;
   REigenVec[2][3] = -Qf * beta_z;
   REigenVec[2][4] =  REigenVec[2][0] * a2;
   REigenVec[2][5] = -Af * beta_y;
   REigenVec[2][6] = -Af * beta_z;

   REigenVec[4][0] =  REigenVec[2][0];
   REigenVec[4][1] = -REigenVec[2][1];
   REigenVec[4][2] = -REigenVec[2][2];
   REigenVec[4][3] = -REigenVec[2][3];
   REigenVec[4][4] =  REigenVec[2][4];
   REigenVec[4][5] =  REigenVec[2][5];
   REigenVec[4][6] =  REigenVec[2][6];

   REigenVec[5][2] = -REigenVec[1][2];
   REigenVec[5][3] = -REigenVec[1][3];
   REigenVec[5][5] =  REigenVec[1][5];
   REigenVec[5][6] =  REigenVec[1][6];

   REigenVec[6][0] =  REigenVec[0][0];
   REigenVec[6][1] = -REigenVec[0][1];
   REigenVec[6][2] = -REigenVec[0][2];
   REigenVec[6][3] = -REigenVec[0][3];
   REigenVec[6][4] =  REigenVec[0][4];
   REigenVec[6][5] =  REigenVec[0][5];
   REigenVec[6][6] =  REigenVec[0][6];


// left eigenvectors (rows of the matrix LEigenVec)
   const real N         = (real)0.5 * _a2;
   const real N_By      = N * beta_y;
   const real N_Bz      = N * beta_z;
   const real As_Rho    = As * _Rho;
   const real Af_Rho    = Af * _Rho;
   const real S_inv_Rho = S * _sqrt_Rho;

   LEigenVec[0][1] = -N * Cff;
   LEigenVec[0][2] =  N_By * Qs;
   LEigenVec[0][3] =  N_Bz * Qs;
   LEigenVec[0][4] =  N * alpha_f * _Rho;
   LEigenVec[0][5] =  N_By * As_Rho;
   LEigenVec[0][6] =  N_Bz * As_Rho;

   LEigenVec[1][2] = -(real)0.5 * beta_z;
   LEigenVec[1][3] =  (real)0.5 * beta_y;
   LEigenVec[1][5] =  LEigenVec[1][2] * S_inv_Rho;
   LEigenVec[1][6] =  LEigenVec[1][3] * S_inv_Rho;

   LEigenVec[2][1] = -N * Css;
   LEigenVec[2][2] = -N_By * Qf;
   LEigenVec[2][3] = -N_Bz * Qf;
   LEigenVec[2][4] =  N * alpha_s * _Rho;
   LEigenVec[2][5] = -N_By * Af_Rho;
   LEigenVec[2][6] = -N_Bz * Af_Rho;

   LEigenVec[3][4] = -_a2;

   LEigenVec[4][1] = -LEigenVec[2][1];
   LEigenVec[4][2] = -LEigenVec[2][2];
   LEigenVec[4][3] = -LEigenVec[2][3];
   LEigenVec[4][4] =  LEigenVec[2][4];
   LEigenVec[4][5] =  LEigenVec[2][5];
   LEigenVec[4][6] =  LEigenVec[2][6];

   LEigenVec[5][2] = -LEigenVec[1][2];
   LEigenVec[5][3] = -LEigenVec[1][3];
   LEigenVec[5][5] =  LEigenVec[1][5];
   LEigenVec[5][6] =  LEigenVec[1][6];

   LEigenVec[6][1] = -LEigenVec[0][1];
   LEigenVec[6][2] = -LEigenVec[0][2];
   LEigenVec[6][3] = -LEigenVec[0][3];
   LEigenVec[6][4] =  LEigenVec[0][4];
   LEigenVec[6][5] =  LEigenVec[0][5];
   LEigenVec[6][6] =  LEigenVec[0][6];


// b. pure hydro
#  else // #ifdef MHD

   const real vx = CC_Var[1];
   const real vy = CC_Var[2];
   const real vz = CC_Var[3];

// eigenvalues along all three spatial directions
   EigenVal[0][0] = vx - a;
   EigenVal[0][1] = vx;
   EigenVal[0][2] = vx;
   EigenVal[0][3] = vx;
   EigenVal[0][4] = vx + a;

   EigenVal[1][0] = vy - a;
   EigenVal[1][1] = vy;
   EigenVal[1][2] = vy;
   EigenVal[1][3] = vy;
   EigenVal[1][4] = vy + a;

   EigenVal[2][0] = vz - a;
   EigenVal[2][1] = vz;
   EigenVal[2][2] = vz;
   EigenVal[2][3] = vz;
   EigenVal[2][4] = vz + a;


// NOTE : the left and right eigenvectors have the same form along different spatial directions for hydro
// left eigenvectors (rows of the matrix LEigenVec)
   LEigenVec[0][1] = -(real)0.5*Rho*_a;
   LEigenVec[0][4] = (real)0.5*_a2;
   LEigenVec[1][4] = -_a2;
   LEigenVec[4][1] = -LEigenVec[0][1];
   LEigenVec[4][4] = +LEigenVec[0][4];


// right eigenvectors (rows instead of columns of the matrix REigenVec for better performance)
   REigenVec[0][1] = -a*_Rho;
   REigenVec[0][4] = a2;
   REigenVec[4][1] = -REigenVec[0][1];
   REigenVec[4][4] = a2;

#  endif // #ifdef MHD ... else ...

} // FUNCTION : Hydro/MHD_GetEigenSystem
#endif // #if (  FLU_SCHEME == CTU  ||  ( defined MHD && defined CHAR_RECONSTRUCTION )  )



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_LimitSlope
// Description :  Evaluate the monotonic slope by slope limiters
//
// Note        :  1. Input data must be primitive variables
//                2. Size of each input array should be NCOMP_LR
//
// Parameter   :  L             : Element x-1
//                C             : Element x
//                R             : Element x+1
//                LR_Limiter    : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff  : Coefficient of the generalized MinMod limiter
//                XYZ           : Target spatial direction : (0/1/2) --> (x/y/z)
//                                --> For CHAR_RECONSTRUCTION only
//                L/REigenVec   : Array storing the left/right eigenvectors
//                                --> For MHD + CHAR_RECONSTRUCTION only
//                Slope_Limiter : Array to store the output monotonic slope
//                EoS           : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_LimitSlope( const real L[], const real C[], const real R[], const LR_Limiter_t LR_Limiter,
                       const real MinMod_Coeff, const int XYZ,
                       const real LEigenVec[][NWAVE], const real REigenVec[][NWAVE], real Slope_Limiter[],
                       const EoS_t *EoS )
{

// check
#  ifdef GAMER_DEBUG
#  if ( defined MHD  &&  defined CHAR_RECONSTRUCTION )
   if ( LEigenVec == NULL )   printf( "ERROR : LEigenVec == NULL !!\n" );
   if ( REigenVec == NULL )   printf( "ERROR : REigenVec == NULL !!\n" );
#  endif
#  endif


   real Slope_L[NCOMP_LR], Slope_R[NCOMP_LR], Slope_C[NCOMP_LR];
   real Slope_A[NCOMP_LR], Slope_LR;

// evaluate different slopes
   for (int v=0; v<NCOMP_LR; v++)
   {
      Slope_L[v] = C[v] - L[v];
      Slope_R[v] = R[v] - C[v];
      Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );
   }

   if ( LR_Limiter == LR_LIMITER_VL_GMINMOD )
   {
      for (int v=0; v<NCOMP_LR; v++)
      {
         if ( Slope_L[v]*Slope_R[v] > (real)0.0 )
            Slope_A[v] = (real)2.0*Slope_L[v]*Slope_R[v]/( Slope_L[v] + Slope_R[v] );
         else
            Slope_A[v] = (real)0.0;
      }
   }


// primitive variables --> characteristic variables
#  ifdef CHAR_RECONSTRUCTION
   const real Dens = C[0];
   const real Pres = C[4];

   Hydro_Pri2Char( Slope_L, Dens, Pres, LEigenVec, XYZ, EoS );
   Hydro_Pri2Char( Slope_R, Dens, Pres, LEigenVec, XYZ, EoS );
   Hydro_Pri2Char( Slope_C, Dens, Pres, LEigenVec, XYZ, EoS );

   if ( LR_Limiter == LR_LIMITER_VL_GMINMOD )
   Hydro_Pri2Char( Slope_A, Dens, Pres, LEigenVec, XYZ, EoS );
#  endif


// apply the slope limiter
   for (int v=0; v<NCOMP_LR; v++)
   {
      Slope_LR = Slope_L[v]*Slope_R[v];

      if ( Slope_LR > (real)0.0 )
      {
         switch ( LR_Limiter )
         {
//          notes for LR_LIMITER_CENTRAL:
//          (1) not TVD --> extra monotonicity check outside this function is required
//          (2) mainly for MHM_RP+PPM to achieve 2nd-order accuracy in linear wave tests
            case LR_LIMITER_CENTRAL:      // central
               Slope_Limiter[v] = Slope_C[v];
               break;

            case LR_LIMITER_VANLEER:      // van-Leer
               Slope_Limiter[v] = (real)2.0*Slope_LR/( Slope_L[v] + Slope_R[v] );
               break;

            case LR_LIMITER_GMINMOD:      // generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            case LR_LIMITER_ALBADA:       // van-Albada
               Slope_Limiter[v] = Slope_LR*( Slope_L[v] + Slope_R[v] ) /
                                  ( Slope_L[v]*Slope_L[v] + Slope_R[v]*Slope_R[v] );
               break;

            case LR_LIMITER_VL_GMINMOD:   // van-Leer + generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_A[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            default :
#              ifdef GAMER_DEBUG
               printf( "ERROR : incorrect parameter %s = %d !!\n", "LR_Limiter", LR_Limiter );
#              endif
               return;
         }
      } // if ( Slope_LR > (real)0.0 )

      else
      {
         Slope_Limiter[v] = (real)0.0;
      } // if ( Slope_LR > (real)0.0 ) ... else ...
   } // for (int v=0; v<NCOMP_LR; v++)


// characteristic variables --> primitive variables
#  ifdef CHAR_RECONSTRUCTION
   Hydro_Char2Pri( Slope_Limiter, Dens, Pres, REigenVec, XYZ, EoS );
#  endif

} // FUNCTION : Hydro_LimitSlope



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_HancockPredict
// Description :  Evolve the face-centered variables by half time-step by calculating the face-centered fluxes
//                (no Riemann solver is required)
//
// Note        :  1. Work for the MHM scheme
//                2. Do NOT require data in the neighboring cells
//                3. Input variables must be conserved variables
//
// Parameter   :  fc                : Face-centered conserved variables to be updated
//                dt                : Time interval to advance solution
//                dh                : Cell size
//                g_cc_array        : Array storing the cell-centered conserved variables for checking
//                                    negative density and pressure
//                                    --> It is just the input array Flu_Array_In[]
//                cc_idx            : Index for accessing g_cc_array[]
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                EoS               : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_HancockPredict( real fc[][NCOMP_LR], const real dt, const real dh,
                           const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                           const real MinDens, const real MinPres, const real MinEint,
                           const EoS_t *EoS )
{

   const real dt_dh2 = (real)0.5*dt/dh;

   real Flux[6][NCOMP_TOTAL], dFlux;


// calculate flux
   for (int f=0; f<6; f++)
      Hydro_Con2Flux( f/2, Flux[f], fc[f], MinPres, EoS->DensEint2Pres_FuncPtr,
                      EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table, NULL );

// update the face-centered variables
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      dFlux = dt_dh2*( Flux[1][v] - Flux[0][v] + Flux[3][v] - Flux[2][v] + Flux[5][v] - Flux[4][v] );

      for (int f=0; f<6; f++)    fc[f][v] -= dFlux;
   }

// check the negative density and energy
   for (int f=0; f<6; f++)
   {
#     ifdef BAROTROPIC_EOS
      if ( fc[f][0] <= (real)0.0 )
#     else
      if ( fc[f][0] <= (real)0.0  ||  fc[f][4] <= (real)0.0 )
#     endif
      {
//       set to the cell-centered values before update
         for (int f=0; f<6; f++)
         for (int v=0; v<NCOMP_TOTAL; v++)
            fc[f][v] = g_cc_array[v][cc_idx];

         break;
      }
   }

// apply density and internal energy floors
   for (int f=0; f<6; f++)
   {
      fc[f][0] = FMAX( fc[f][0], MinDens );
#     ifndef BAROTROPIC_EOS
#     ifdef MHD
#     error : ERROR : MHD is not supported here !!!
      const real Emag = NULL_REAL;
#     else
      const real Emag = NULL_REAL;
#     endif
      fc[f][4] = Hydro_CheckMinEintInEngy( fc[f][0], fc[f][1], fc[f][2], fc[f][3], fc[f][4],
                                           MinEint, Emag );
#     endif // #ifndef BAROTROPIC_EOS
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      fc[f][v] = FMAX( fc[f][v], TINY_NUMBER );
#     endif
   }

} // FUNCTION : Hydro_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



// MINMOD macro is only used in this function
#ifdef MINMOD
#  undef MINMOD
#endif



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



#endif // #ifndef __CUFLU_DATARECONSTRUCTION__
