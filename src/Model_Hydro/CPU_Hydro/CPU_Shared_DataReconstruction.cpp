#include "GAMER.h"
#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



extern void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward );
extern real CPU_CheckMinPres( const real InPres, const real MinPres );

static void Get_EigenSystem( const real CC_Var[], real EigenVal[][NCOMP_FLUID], real LEigenVec[][NCOMP_FLUID],
                             real REigenVec[][NCOMP_FLUID], const real Gamma );
static void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                        const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                        const real Gamma, const int XYZ, real Slope_Limiter[] );
extern void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                         const bool NormPassive, const int NNorm, const int NormIdx[] );
#if ( FLU_SCHEME == MHM )
extern void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres );
static void CPU_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                                const real Gamma_m1, const real _Gamma_m1,
                                const real cc_array[][ FLU_NXT*FLU_NXT*FLU_NXT ], const int cc_idx,
                                const real MinDens, const real MinPres );
#endif
#ifdef CHAR_RECONSTRUCTION
static void Pri2Char( real Var[], const real Gamma, const real Rho, const real Pres, const int XYZ );
static void Char2Pri( real Var[], const real Gamma, const real Rho, const real Pres, const int XYZ );
#endif




#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. Input data are primitive variables while output data are conserved variables
//                   --> The latter is because CPU_HancockPredict() works with conserved variables
//                3. PLM and PPM data reconstruction functions share the same function name
//                4. Face-centered variables will be advanced by half time-step for the MHM and CTU schemes
//                5. Data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION" in the header CUFLU.h
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//
//
// Parameter   :  PriVar         : Array storing the input cell-centered primitive variables
//                ConVar         : Array storing the input cell-centered conserved variables
//                                 --> Only for checking negative density and pressure in CPU_HancockPredict()
//                FC_Var         : Array to store the output face-centered primitive variables
//                NIn            : Size of PriVar[] along each direction
//                                 --> Can be smaller than FLU_NXT
//                NGhost         : Number of ghost zones
//                                  --> "NIn-2*NGhost" cells will be computed along each direction
//                                  --> Size of FC_Var[] is assumed to be "(NIn-2*NGhost)^3"
//                                  --> The reconstructed data at cell (i,j,k) will be stored in the
//                                      array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma          : Ratio of specific heats
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                dt             : Time interval to advance solution (for the CTU scheme)
//                dh             : Cell size
//                MinDens/Pres   : Minimum allowed density and pressure
//                NormPassive    : true --> convert passive scalars to mass fraction
//                NNorm          : Number of passive scalars for the option "NormPassive"
//                                 --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx        : Target variable indices for the option "NormPassive"
//                                 --> Should be set to the global variable "PassiveNorm_VarIdx"
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[][ FLU_NXT*FLU_NXT*FLU_NXT    ],
                             const real ConVar[][ FLU_NXT*FLU_NXT*FLU_NXT    ],
                                   real FC_Var[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                             const int NIn, const int NGhost, const real Gamma,
                             const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh,
                             const real MinDens, const real MinPres,
                             const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const int  didx_cc[3] = { 1, NIn, NIn*NIn };
   const int  NOut       = NIn - 2*NGhost;      // number of output cells
   const real  Gamma_m1  = Gamma - (real)1.0;
   const real _Gamma_m1  = (real)1.0 / Gamma_m1;

// variables for the CTU scheme
   /*
#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt/dh;

   real EigenVal[3][NCOMP_FLUID], Correct_L[NCOMP_TOTAL], Correct_R[NCOMP_TOTAL], dFC[NCOMP_TOTAL];
   real Coeff_L, Coeff_R;

// initialize the constant components of the matrices of the left and right eigenvectors
   real LEigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #if ( FLU_SCHEME ==  CTU )
   */


#  ifdef __CUDACC__
   for (int idx_fc=threadIdx.x; idx_fc<CUBE(NOut); idx_fc+=blockDim.x)
   {
      const int size_ij = SQR(NOut);
      const int i_cc    = NGhost + idx_fc%NOut;
      const int j_cc    = NGhost + idx_fc%size_ij/NOut;
      const int k_cc    = NGhost + idx_fc/size_ij;
#  else
   for (int k_cc=NGhost, k_fc=0;  k_cc<NGhost+NOut;  k_cc++, k_fc++)
   for (int j_cc=NGhost, j_fc=0;  j_cc<NGhost+NOut;  j_cc++, j_fc++)
   for (int i_cc=NGhost, i_fc=0;  i_cc<NGhost+NOut;  i_cc++, i_fc++)
   {
      const int idx_cc = IDX321( i_cc, j_cc, k_cc, NIn,  NIn  );
      const int idx_fc = IDX321( i_fc, j_fc, k_fc, NOut, NOut );
#  endif

      real cc_C[NCOMP_TOTAL], cc_L[NCOMP_TOTAL], cc_R[NCOMP_TOTAL];  // cell-centered variables of the Central/Left/Right cells
      real fc[6][NCOMP_TOTAL];                                       // face-centered variables of the central cell
      real Slope_Limiter[NCOMP_TOTAL];

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C[v] = PriVar[v][idx_cc];


//    1. evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Get_EigenSystem( cc_C, EigenVal, LEigenVec, REigenVec, Gamma );
#     endif


//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
//       2. evaluate the monotonic slope
         const int faceL   = 2*d;      // left and right face indices
         const int faceR   = faceL+1;
         const int idx_ccL = idx_cc - didx_cc[d];
         const int idx_ccR = idx_cc + didx_cc[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            cc_L[v] = PriVar[v][idx_ccL];
            cc_R[v] = PriVar[v][idx_ccR];
         }

         LimitSlope( NULL, cc_L, cc_C, cc_R, NULL, LR_Limiter, MinMod_Coeff, NULL_REAL, Gamma, d, Slope_Limiter );


//       3. get the face-centered primitive variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fc[faceL][v] = cc_C[v] - (real)0.5*Slope_Limiter[v];
            fc[faceR][v] = cc_C[v] + (real)0.5*Slope_Limiter[v];
         }

//       ensure the face-centered variables lie between neighboring cell-centered values
         for (int v=0; v<NCOMP_TOTAL; v++)
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
         /*
#        if ( FLU_SCHEME == CTU )

//       4-1. evaluate the slope (for passive scalars as well)
         for (int v=0; v<NCOMP_TOTAL; v++)   dFC[v] = fc[faceR][v] - fc[faceL][v];


//       4-2. re-order variables for the y/z directions
         CPU_Rotate3D( dFC, d, true );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

//       4-2-a1. evaluate the corrections to the left and right face-centered variables

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = (real)0.0;
            Correct_R[v] = (real)0.0;
         }

#        ifdef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;

            for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dFC[v];

            Coeff_L *= -dt_dh2*EigenVal[d][Mode];

            for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

         for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] = Correct_L[v];

#        else // ifndef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dFC[v];

               Coeff_L *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*dFC[v];

               Coeff_R *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // ifdef HLL_INCLUDE_ALL_WAVES ... else ...


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT || ifndef HLL_NO_REF_STATE )

//       4-2-b1. evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][0], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][4], (real)0.0 );

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = Coeff_L*dFC[v];
            Correct_R[v] = Coeff_R*dFC[v];
         }


//       4-2-b2. evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dFC[v];

               Coeff_L *= dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*dFC[v];

               Coeff_R *= dt_dh2*( EigenVal[d][4] - EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       4-2-b3. evaluate the corrections to the left and right face-centered passive scalars
//               --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*dFC[v];
            Correct_R[v] = Coeff_R*dFC[v];
         }
#        endif


//       4-2-b4. evaluate the face-centered variables at the half time-step
         CPU_Rotate3D( Correct_L, d, false );
         CPU_Rotate3D( Correct_R, d, false );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fc[faceL][v] += Correct_L[v];
            fc[faceR][v] += Correct_R[v];
         }

//       ensure positive density and pressure
         fc[faceL][0] = FMAX( fc[faceL][0], MinDens );
         fc[faceR][0] = FMAX( fc[faceR][0], MinDens );

         fc[faceL][4] = CPU_CheckMinPres( fc[faceL][4], MinPres );
         fc[faceR][4] = CPU_CheckMinPres( fc[faceR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         fc[faceL][v] = FMAX( fc[faceL][v], TINY_NUMBER );
         fc[faceR][v] = FMAX( fc[faceR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )
         */


//       5. primitive variables --> conserved variables
         real tmp[NCOMP_TOTAL];  // input and output arrays must not overlap for Pri2Con()

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceL][v];
         CPU_Pri2Con( tmp, fc[faceL], _Gamma_m1, NormPassive, NNorm, NormIdx );

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceR][v];
         CPU_Pri2Con( tmp, fc[faceR], _Gamma_m1, NormPassive, NNorm, NormIdx );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    6. advance the face-centered variables by half time-step for the MHM integrator
//       --> here we have assumed that ConVar[]
      CPU_HancockPredict( fc, dt, dh, Gamma_m1, _Gamma_m1, ConVar, idx_cc, MinDens, MinPres );
#     endif


//    7. store the face-centered values to the output array
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         FC_Var[f][v][idx_fc] = fc[f][v];

   } // k,j,i

} // FUNCTION : CPU_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. Input data are primitive variables while output data are conserved variables
//                   --> The latter is because CPU_HancockPredict() works with conserved variables
//                3. PLM and PPM data reconstruction functions share the same function name
//                4. Face-centered variables will be advanced by half time-step for the MHM and CTU schemes
//                5. Data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION" in the header CUFLU.h
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  PriVar         : Array storing the input cell-centered primitive variables
//                FC_Var         : Array to store the output face-centered primitive variables
//                NIn            : Size of PriVar[] along each direction
//                NGhost         : Number of ghost zones
//                                  --> "NIn-2*NGhost" cells will be computed along each direction
//                                  --> Size of FC_Var[] is assumed to be "(NIn-2*NGhost)^3"
//                                  --> The reconstructed data at cell (i,j,k) will be stored in the
//                                      array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma          : Ratio of specific heats
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                dt             : Time interval to advance solution (for the CTU scheme)
//                dh             : Cell size (for the CTU scheme)
//                MinDens/Pres   : Minimum allowed density and pressure
//                NormPassive    : true --> convert passive scalars to mass fraction
//                NNorm          : Number of passive scalars for the option "NormPassive"
//                                 --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx        : Target variable indices for the option "NormPassive"
//                                 --> Should be set to the global variable "PassiveNorm_VarIdx"
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[]   [ FLU_NXT*FLU_NXT*FLU_NXT    ],
                                   real FC_Var[][6][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                             const int NIn, const int NGhost, const real Gamma,
                             const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh,
                             const real MinDens, const real MinPres,
                             const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const int NOut   = NIn - 2*NGhost;                    // number of output cells
   const int NSlope = NOut + 2;                          // number of cells required to store the slope data
   const int didx_cc[3] = { 1, NIn, NIn*NIn };
   const int dr3[3] = { 1, NSlope, NSlope*NSlope };

   int ID1, ID2, ID3, ID1_L, ID1_R, ID3_L, ID3_R, dL, dR;
   real Slope_Limiter[NCOMP_TOTAL] = { (real)0.0 };
   real CC_L, CC_R, CC_C, dCC_L, dCC_R, dCC_C, FC_L, FC_R, dFC[NCOMP_TOTAL], dFC6[NCOMP_TOTAL], Max, Min;

   real (*Slope_PPM)[3][NCOMP_TOTAL] = new real [ NSlope*NSlope*NSlope ][3][NCOMP_TOTAL];

// variables for the CTU scheme
#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt/dh;

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )
#  ifdef HLL_INCLUDE_ALL_WAVES
   const bool HLL_Include_All_Waves = true;
#  else
   const bool HLL_Include_All_Waves = false;
#  endif
#  endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

   real EigenVal[3][NCOMP_FLUID], Correct_L[NCOMP_TOTAL], Correct_R[NCOMP_TOTAL];
   real Coeff_A, Coeff_B, Coeff_C, Coeff_D, Coeff_L, Coeff_R;

// initialize the constant components of the matrices of the left and right eigenvectors
   real LEigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #if ( FLU_SCHEME == CTU )


// (2-1) evaluate the monotonic slope
   for (int k1=NGhost-1, k2=0;  k1<NGhost-1+NSlope;  k1++, k2++)
   for (int j1=NGhost-1, j2=0;  j1<NGhost-1+NSlope;  j1++, j2++)
   for (int i1=NGhost-1, i2=0;  i1<NGhost-1+NSlope;  i1++, i2++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NSlope + j2)*NSlope + i2;

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         ID1_L = ID1 - didx_cc[d];
         ID1_R = ID1 + didx_cc[d];

         LimitSlope( NULL, PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], NULL, LR_Limiter,
                     MinMod_Coeff, NULL_REAL, Gamma, d, Slope_Limiter );

//       store the slope to the array "Slope_PPM"
         for (int v=0; v<NCOMP_TOTAL; v++)   Slope_PPM[ID2][d][v] = Slope_Limiter[v];

      } // for (int d=0; d<3; d++)
   } // k,j,i


   for (int k1=NGhost, k2=0, k3=1;  k1<NGhost+NOut;  k1++, k2++, k3++)
   for (int j1=NGhost, j2=0, j3=1;  j1<NGhost+NOut;  j1++, j2++, j3++)
   for (int i1=NGhost, i2=0, i3=1;  i1<NGhost+NOut;  i1++, i2++, i3++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NOut   + j2)*NOut   + i2;
      ID3 = (k3*NSlope + j3)*NSlope + i3;


//    (2-2) evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Get_EigenSystem( PriVar[ID1], EigenVal, LEigenVec, REigenVec, Gamma );
#     endif


//    (2-3) get the face-centered primitive variables
//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         dL    = 2*d;
         dR    = dL+1;
         ID1_L = ID1 - didx_cc[d];
         ID1_R = ID1 + didx_cc[d];
         ID3_L = ID3 - dr3[d];
         ID3_R = ID3 + dr3[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          (2-3-1) parabolic interpolation
            CC_L  = PriVar[ID1_L][v];
            CC_R  = PriVar[ID1_R][v];
            CC_C  = PriVar[ID1  ][v];

            dCC_L = Slope_PPM[ID3_L][d][v];
            dCC_R = Slope_PPM[ID3_R][d][v];
            dCC_C = Slope_PPM[ID3  ][d][v];

            FC_L  = (real)0.5*( CC_C + CC_L ) - (real)1.0/(real)6.0*( dCC_C - dCC_L );
            FC_R  = (real)0.5*( CC_C + CC_R ) - (real)1.0/(real)6.0*( dCC_R - dCC_C );


//          (2-3-2) monotonicity constraint
            dFC [v] = FC_R - FC_L;
            dFC6[v] = (real)6.0*(  CC_C - (real)0.5*( FC_L + FC_R )  );

            if (  ( FC_R - CC_C )*( CC_C - FC_L ) <= (real)0.0  )
            {
               FC_L = CC_C;
               FC_R = CC_C;
            }
            else if ( dFC[v]*dFC6[v] > +dFC[v]*dFC[v] )
               FC_L = (real)3.0*CC_C - (real)2.0*FC_R;
            else if ( dFC[v]*dFC6[v] < -dFC[v]*dFC[v] )
               FC_R = (real)3.0*CC_C - (real)2.0*FC_L;


//          (2-3-3) ensure the face-centered variables lie between neighboring cell-centered values
            Min  = ( CC_C < CC_L ) ? CC_C : CC_L;
            Max  = ( CC_C > CC_L ) ? CC_C : CC_L;
            FC_L = ( FC_L > Min  ) ? FC_L : Min;
            FC_L = ( FC_L < Max  ) ? FC_L : Max;

            Min  = ( CC_C < CC_R ) ? CC_C : CC_R;
            Max  = ( CC_C > CC_R ) ? CC_C : CC_R;
            FC_R = ( FC_R > Min  ) ? FC_R : Min;
            FC_R = ( FC_R < Max  ) ? FC_R : Max;


            FC_Var[ID2][dL][v] = FC_L;
            FC_Var[ID2][dR][v] = FC_R;

         } // for (int v=0; v<NCOMP_TOTAL; v++)


//       (2-4) advance the face-centered variables by half time-step for the CTU integrator
#        if ( FLU_SCHEME == CTU )

//       (2-4-1) compute the PPM coefficient (for the passive scalars as well)
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            dFC [v] = FC_Var[ID2][dR][v] - FC_Var[ID2][dL][v];
            dFC6[v] = (real)6.0*(  PriVar[ID1][v] - (real)0.5*( FC_Var[ID2][dL][v] + FC_Var[ID2][dR][v] )  );
         }

//       (2-4-2) re-order variables for the y/z directions
         CPU_Rotate3D( dFC,  d, true );
         CPU_Rotate3D( dFC6, d, true );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

//       (2-4-a1) evaluate the corrections to the left and right face-centered variables
         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = (real)0.0;
            Correct_R[v] = (real)0.0;
         }

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] <= (real)0.0 )
            {
               Coeff_C = -dt_dh2*EigenVal[d][Mode];
               Coeff_D = real(-4.0/3.0)*SQR(Coeff_C);

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dFC[v] + dFC6[v] ) +
                                                                                    Coeff_D*( dFC6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] >= (real)0.0 )
            {
               Coeff_A = -dt_dh2*EigenVal[d][Mode];
               Coeff_B = real(-4.0/3.0)*SQR(Coeff_A);

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dFC[v] - dFC6[v] ) +
                                                                                    Coeff_B*( dFC6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT && ifndef HLL_NO_REF_STATE )

//       (2-4-b1) evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][0], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][4], (real)0.0 );

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = Coeff_L*(  dFC[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6[v]  );
            Correct_R[v] = Coeff_R*(  dFC[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6[v]  );
         }


//       (2-4-b2) evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               Coeff_C = dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_D=0 when Coeff_C=0
/*             Coeff_D = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][   0]*EigenVal[d][   0] -
                                                        EigenVal[d][Mode]*EigenVal[d][Mode]   ); */
               Coeff_D = real(4.0/3.0)*dt_dh2*Coeff_C*( EigenVal[d][0] + EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dFC[v] + dFC6[v] ) +
                                                                                    Coeff_D*( dFC6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               Coeff_A = dt_dh2*( EigenVal[d][4] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_B=0 when Coeff_A=0
/*             Coeff_B = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][   4]*EigenVal[d][   4] -
                                                        EigenVal[d][Mode]*EigenVal[d][Mode]   ); */
               Coeff_B = real(4.0/3.0)*dt_dh2*Coeff_A*( EigenVal[d][4] + EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dFC[v] - dFC6[v] ) +
                                                                                    Coeff_B*( dFC6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       (2-4-3) evaluate the corrections to the left and right face-centered passive scalars
//               --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*(  dFC[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6[v]  );
            Correct_R[v] = Coeff_R*(  dFC[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6[v]  );
         }
#        endif


//       (2-4-4) evaluate the face-centered variables at the half time-step
         CPU_Rotate3D( Correct_L, d, false );
         CPU_Rotate3D( Correct_R, d, false );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            FC_Var[ID2][dL][v] += Correct_L[v];
            FC_Var[ID2][dR][v] += Correct_R[v];
         }

//       ensure positive density and pressure
         FC_Var[ID2][dL][0] = FMAX( FC_Var[ID2][dL][0], MinDens );
         FC_Var[ID2][dR][0] = FMAX( FC_Var[ID2][dR][0], MinDens );

         FC_Var[ID2][dL][4] = CPU_CheckMinPres( FC_Var[ID2][dL][4], MinPres );
         FC_Var[ID2][dR][4] = CPU_CheckMinPres( FC_Var[ID2][dR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         FC_Var[ID2][dL][v] = FMAX( FC_Var[ID2][dL][v], TINY_NUMBER );
         FC_Var[ID2][dR][v] = FMAX( FC_Var[ID2][dR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )

      } // for (int d=0; d<3; d++)
   } // k,j,i

   delete [] Slope_PPM;


   CONSERVE --> PREMITIVE


} // FUNCTION : CPU_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )



/*
#ifdef CHAR_RECONSTRUCTION
//-------------------------------------------------------------------------------------------------------
// Function    :  Pri2Char
// Description :  Convert the primitive variables to the characteristic variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  InOut : Array storing both the input primitive variables and output characteristic variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
void Pri2Char( real InOut[], const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(Pres) )
      Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( CPU_CheckNegative(Rho) )
      Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Cs2 = (real)1.0 / ( Gamma*Pres/Rho );
   const real _Cs  = SQRT( _Cs2 );
   real Temp[NCOMP_FLUID];

   for (int v=0; v<NCOMP_FLUID; v++)   Temp[v] = InOut[v];

   CPU_Rotate3D( Temp, XYZ, true );

   InOut[0] = -(real)0.5*Rho*_Cs*Temp[1] + (real)0.5*_Cs2*Temp[4];
   InOut[1] = Temp[0] - _Cs2*Temp[4];
   InOut[2] = Temp[2];
   InOut[3] = Temp[3];
   InOut[4] = +(real)0.5*Rho*_Cs*Temp[1] + (real)0.5*_Cs2*Temp[4];

} // FUNCTION : Pri2Char



//-------------------------------------------------------------------------------------------------------
// Function    :  Char2Pri
// Description :  Convert the characteristic variables to the primitive variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  InOut : Array storing both the input characteristic variables and output primitive variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
void Char2Pri( real InOut[], const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(Pres) )
      Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( CPU_CheckNegative(Rho) )
      Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Rho = (real)1.0 / Rho;
   const real Cs2  = Gamma*Pres*_Rho;
   const real Cs   = SQRT( Cs2 );
   real Temp[NCOMP_FLUID];

   Temp[0] = InOut[0] + InOut[1] + InOut[4];
   Temp[1] = Cs*_Rho*( -InOut[0] + InOut[4] );
   Temp[2] = InOut[2];
   Temp[3] = InOut[3];
   Temp[4] = Cs2*( InOut[0] + InOut[4] );

   CPU_Rotate3D( Temp, XYZ, false );

   for (int v=0; v<NCOMP_FLUID; v++)   InOut[v] = Temp[v];

} // FUNCTION : Char2Pri
#endif
*/



/*
#if ( FLU_SCHEME == CTU )
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_EigenSystem
// Description :  Evaluate the eigenvalues and left/right eigenvectors
//
// Note        :  1. The input data should be primitive variables
//                2. The constant components of eigenvectors should be initialized in advance
//                3. Work for the CTU scheme
//                4. Do not need to consider passive scalars
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  CC_Var      : Array storing the input cell-centered primitive variables
//                EigenVal    : Array to store the output eigenvalues (in three spatial directions)
//                LEigenVec   : Array to store the output left eigenvectors
//                REigenVec   : Array to store the output right eigenvectors
//                Gamma       : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void Get_EigenSystem( const real CC_Var[], real EigenVal[][NCOMP_FLUID], real LEigenVec[][NCOMP_FLUID],
                      real REigenVec[][NCOMP_FLUID], const real Gamma )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(CC_Var[4]) )
      Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   CC_Var[4], __FILE__, __LINE__, __FUNCTION__ );

   if ( CPU_CheckNegative(CC_Var[0]) )
      Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   CC_Var[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real  Rho = CC_Var[0];
   const real _Rho = (real)1.0/Rho;
   const real  Vx  = CC_Var[1];
   const real  Vy  = CC_Var[2];
   const real  Vz  = CC_Var[3];
   const real Pres = CC_Var[4];
   const real  Cs2 = Gamma*Pres*_Rho;
   const real  Cs  = SQRT( Cs2 );
   const real _Cs  = (real)1.0/Cs;
   const real _Cs2 = _Cs*_Cs;

// a. eigenvalues in three spatial directions
   EigenVal[0][0] = Vx - Cs;
   EigenVal[0][1] = Vx;
   EigenVal[0][2] = Vx;
   EigenVal[0][3] = Vx;
   EigenVal[0][4] = Vx + Cs;

   EigenVal[1][0] = Vy - Cs;
   EigenVal[1][1] = Vy;
   EigenVal[1][2] = Vy;
   EigenVal[1][3] = Vy;
   EigenVal[1][4] = Vy + Cs;

   EigenVal[2][0] = Vz - Cs;
   EigenVal[2][1] = Vz;
   EigenVal[2][2] = Vz;
   EigenVal[2][3] = Vz;
   EigenVal[2][4] = Vz + Cs;


// NOTE : the left and right eigenvectors have the same form in different spatial directions
// b. left eigenvectors (rows of the matrix LEigenVec)
   LEigenVec[0][1] = -(real)0.5*Rho*_Cs;
   LEigenVec[0][4] = (real)0.5*_Cs2;
   LEigenVec[1][4] = -_Cs2;
   LEigenVec[4][1] = -LEigenVec[0][1];
   LEigenVec[4][4] = +LEigenVec[0][4];


// c. right eigenvectors (rows of the matrix REigenVec)
   REigenVec[0][1] = -Cs*_Rho;
   REigenVec[0][4] = Cs2;
   REigenVec[4][1] = -REigenVec[0][1];
   REigenVec[4][4] = Cs2;

} // FUNCTION : Get_EigenSystem
#endif // #if ( FLU_SCHEME == CTU )
*/



//-------------------------------------------------------------------------------------------------------
// Function    :  LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//                2. The L2 and R2 elements are useful only for the extrema-preserving limiter
//
// Parameter   :  L2             : Element x-2
//                L1             : Element x-1
//                C0             : Element x
//                R1             : Element x+1
//                R2             : Element x+2
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                Gamma          : Ratio of specific heats
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ            : Target spatial direction : (0/1/2) --> (x/y/z)
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                Slope_Limiter  : Array to store the output monotonic slope
//-------------------------------------------------------------------------------------------------------
void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                 const real Gamma, const int XYZ, real Slope_Limiter[] )
{

#  ifdef CHAR_RECONSTRUCTION
   const real Rho  = C0[0];
   const real Pres = C0[4];
#  endif

   real Slope_L[NCOMP_TOTAL], Slope_R[NCOMP_TOTAL], Slope_C[NCOMP_TOTAL], Slope_A[NCOMP_TOTAL], Slope_LR;


// evaluate different slopes
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_L[v] = C0[v] - L1[v];
      Slope_R[v] = R1[v] - C0[v];
      Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );
   }

   if ( LR_Limiter == VL_GMINMOD )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         if ( Slope_L[v]*Slope_R[v] > (real)0.0 )
            Slope_A[v] = (real)2.0*Slope_L[v]*Slope_R[v]/( Slope_L[v] + Slope_R[v] );
         else
            Slope_A[v] = (real)0.0;
      }
   }


// primitive variables --> characteristic variables
#  ifdef CHAR_RECONSTRUCTION
   Pri2Char( Slope_L, Gamma, Rho, Pres, XYZ );
   Pri2Char( Slope_R, Gamma, Rho, Pres, XYZ );
   Pri2Char( Slope_C, Gamma, Rho, Pres, XYZ );

   if ( LR_Limiter == VL_GMINMOD )
      Pri2Char( Slope_A, Gamma, Rho, Pres, XYZ );
#  endif


// apply the slope limiter
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_LR = Slope_L[v]*Slope_R[v];

      if ( Slope_LR > (real)0.0 )
      {
         switch ( LR_Limiter )
         {
            case VANLEER:     // van-Leer
               Slope_Limiter[v] = (real)2.0*Slope_LR/( Slope_L[v] + Slope_R[v] );
               break;

            case GMINMOD:     // generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            case ALBADA:      // van-Albada
               Slope_Limiter[v] = Slope_LR*( Slope_L[v] + Slope_R[v] ) /
                                  ( Slope_L[v]*Slope_L[v] + Slope_R[v]*Slope_R[v] );
               break;

            case VL_GMINMOD:  // van-Leer + generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_A[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            default :
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LR_Limiter", LR_Limiter );
         }
      } // if ( Slope_LR > (real)0.0 )

      else
      {
         Slope_Limiter[v] = (real)0.0;
      } // if ( Slope_LR > (real)0.0 ) ... else ...
   } // for (int v=0; v<NCOMP_TOTAL; v++)


// characteristic variables --> primitive variables
#  ifdef CHAR_RECONSTRUCTION
   Char2Pri( Slope_Limiter, Gamma, Rho, Pres, XYZ );
#  endif

} // FUNCTION : LimitSlope



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_HancockPredict
// Description :  Evolve the face-centered variables by half time-step by calculating the face-centered fluxes
//                (no Riemann solver is required)
//
// Note        :  1. Work for the MHM scheme
//                2. Do NOT require data in the neighboring cells
//                3. Input variables must be conserved variables
//
// Parameter   :  fc           : Face-centered conserved variables to be updated
//                dt           : Time interval to advance solution
//                dh           : Cell size
//                Gamma_m1     : Gamma - 1
//                _Gamma_m1    : 1 / (Gamma - 1)
//                cc_array     : Array storing the cell-centered conserved variables for checking
//                               negative density and pressure
//                               --> It is just the input array Flu_Array_In[]
//                cc_idx       : Index for accessing cc_array[]
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                         const real Gamma_m1, const real _Gamma_m1,
                         const real cc_array[][ FLU_NXT*FLU_NXT*FLU_NXT ], const int cc_idx,
                         const real MinDens, const real MinPres )
{

   const real dt_dh2 = (real)0.5*dt/dh;

   real Flux[6][NCOMP_TOTAL], dFlux[NCOMP_TOTAL];


// calculate flux
   for (int f=0; f<6; f++)    CPU_Con2Flux( f/2, Flux[f], fc[f], Gamma_m1, MinPres );

// update the face-centered variables
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      dFlux[v] = dt_dh2*( Flux[1][v] - Flux[0][v] + Flux[3][v] - Flux[2][v] + Flux[5][v] - Flux[4][v] );

      for (int f=0; f<6; f++)    fc[f][v] -= dFlux[v];
   }

// check the negative density and energy
   for (int f=0; f<6; f++)
   {
      if ( fc[f][0] <= (real)0.0  ||  fc[f][4] <= (real)0.0 )
      {
//       set to the cell-centered values before update
         for (int f=0; f<6; f++)
         for (int v=0; v<NCOMP_TOTAL; v++)
            fc[f][v] = cc_array[v][cc_idx];

         break;
      }
   }

// ensure positive density and pressure
   for (int f=0; f<6; f++)
   {
      fc[f][0] = FMAX( fc[f][0], MinDens );
      fc[f][4] = CPU_CheckMinPresInEngy( fc[f][0], fc[f][1], fc[f][2], fc[f][3], fc[f][4],
                                         Gamma_m1, _Gamma_m1, MinPres );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      fc[f][v] = FMAX( fc[f][v], TINY_NUMBER );
#     endif
   }

} // FUNCTION : CPU_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP || CTU) )
