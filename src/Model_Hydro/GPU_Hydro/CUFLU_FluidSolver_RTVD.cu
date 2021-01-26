#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )


// check before compiling anything else
#if ( NCOMP_PASSIVE != 0 )
#  error : RTVD scheme does NOT support passive scalars !!
#endif


#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

#define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
#define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

static __device__ void CUFLU_Advance( real g_Fluid_In [][5][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][5][ CUBE(PS2) ],
                                      real g_Flux[][9][5][ SQR(PS2) ],
                                      const real dt, const real _dh, const bool StoreFlux,
                                      const int j_gap, const int k_gap, real s_cu[][5][FLU_NXT],
                                      real s_cw[][5][FLU_NXT], real s_flux[][5][FLU_NXT], real s_RLflux[][5][FLU_NXT],
                                      const bool FinalOut, const int XYZ,
                                      const real MinDens, const real MinPres, const real MinEint,
                                      const EoS_DE2P_t EoS_DensEint2Pres,
                                      const EoS_DP2C_t EoS_DensPres2CSqr,
                                      const double EoS_AuxArray_Flt[],
                                      const int    EoS_AuxArray_Int[],
                                      const real *const EoS_Table[EOS_NTABLE_MAX] );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_RTVD
// Description :  GPU fluid solver based on the relaxing TVD (RTVD) scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The three-dimensional evolution is achieved by using the dimensional-split method
//
// Parameter   :  g_Fluid_In             : Global memory array to store the input fluid variables
//                g_Fluid_Out            : Global memory array to store the output fluid variables
//                g_Flux                 : Global memory array to store the output fluxes
//                g_Corner               : Global memory array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                g_Pot_USG              : Global memory array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                dt                     : Time interval to advance solution
//                _dh                    : 1 / grid size
//                StoreFlux              : true --> store the coarse-fine fluxes
//                XYZ                    : true  : x->y->z ( forward sweep)
//                                         false : z->y->x (backward sweep)
//                MinDens                : Density floor
//                MinPres                : Pressure floor
//                MinEint                : Internal energy floor
//                EoS_DensEint2Pres_Func : Function pointer to the EoS routine of computing the gas pressure
//                EoS_DensPres2Eint_Func :                    . . .                             gas internal energy
//                EoS_DensPres2CSqr_Func :                    . . .                             sound speed square
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_RTVD(
   real g_Fluid_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real g_Fluid_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real g_Flux     [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner[][3],
   const real g_Pot_USG[][ CUBE(USG_NXT_F) ],
   const real dt, const real _dh, const bool StoreFlux,
   const bool XYZ, const real MinDens, const real MinPres, const real MinEint,
   const EoS_DE2P_t EoS_DensEint2Pres_Func,
   const EoS_DP2E_t EoS_DensPres2Eint_Func,
   const EoS_DP2C_t EoS_DensPres2CSqr_Func )
{

   __shared__ real s_cu    [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_cw    [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_flux  [FLU_BLOCK_SIZE_Y][5][FLU_NXT];
   __shared__ real s_RLflux[FLU_BLOCK_SIZE_Y][5][FLU_NXT];

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux,              0,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 0, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux, FLU_GHOST_SIZE,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 3, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux,  true, 6, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux,              0,              0,
                     s_cu, s_cw, s_flux, s_RLflux, false, 6, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux,              0, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux, false, 3, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );

      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, StoreFlux, FLU_GHOST_SIZE, FLU_GHOST_SIZE,
                     s_cu, s_cw, s_flux, s_RLflux,  true, 0, MinDens, MinPres, MinEint,
                     EoS_DensEint2Pres_Func, EoS_DensPres2CSqr_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );
   }

} // FUNCTION : CUFLU_FluidSolver_RTVD



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  GPU device function, which performs a one-dimensional sweep based on the TVD scheme
//
// Note        :  a. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                b. The direction of the one dimensional sweep is determined by the input parameter "XYZ"
//
// Parameter   :  g_Fluid_In        : Global memory array to store the input fluid variables
//                g_Fluid_Out       : Global memory array to store the output fluid variables
//                g_Flux            : Global memory array to store the output fluxes
//                dt                : Time interval to advance solution
//                _dh               : 1 / grid size
//                StoreFlux         : true --> store the coarse-fine fluxes
//                j_gap             : Number of useless grids in each side in the j direction (j may not be equal to y)
//                k_gap             : Number of useless grids in each side in the k direction (k mya not be equal to z)
//                s_cu              : Shared memory array storing the normal flux
//                s_cw              : Shared memory array storing the auxiliary flux
//                s_flux            : Shared memory array storing the final flux used to update the fluid variables
//                s_RLflux          : Shared memory array storing the left/right-moving flux
//                XYZ               : 0 : Update the solution in the x direction
//                                    3 : Update the solution in the y direction
//                                    6 : Update the solution in the z direction
//                                    --> This parameter is also used to determine the place to store the output fluxes
//                MinDens           : Density floor
//                MinPres           : Pressure floor
//                MinEint           : Internal energy floor
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Advance( real g_Fluid_In [][5][ CUBE(FLU_NXT) ],
                               real g_Fluid_Out[][5][ CUBE(PS2) ],
                               real g_Flux[][9][5][ SQR(PS2) ],
                               const real dt, const real _dh, const bool StoreFlux,
                               const int j_gap, const int k_gap, real s_cu[][5][FLU_NXT],
                               real s_cw[][5][FLU_NXT], real s_flux[][5][FLU_NXT], real s_RLflux[][5][FLU_NXT],
                               const bool FinalOut, const int XYZ,
                               const real MinDens, const real MinPres, const real MinEint,
                               const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr,
                               const double EoS_AuxArray_Flt[],
                               const int    EoS_AuxArray_Int[],
                               const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   const uint bx               = blockIdx.x;
   const uint tx               = threadIdx.x;
   const uint ty               = threadIdx.y;
   const uint dj               = blockDim.y;
   const uint size_j           = FLU_NXT - (j_gap<<1);
   const uint size_k           = FLU_NXT - (k_gap<<1);
   const uint NColumn          = __umul24( size_j, size_k );
   const uint i                = tx;                  // (i,j) the element in shared memory under evaluation
   const uint ip               = i+1;
   const uint im               = i-1;
         uint j                = j_gap + ty%size_j;
         uint k                = k_gap + ty/size_j;
         uint Column0          = 0;                   // the total number of columns that have been updated
   const uint j_end            = FLU_NXT - j_gap;
   const uint k_end            = FLU_NXT - k_gap;
   const real dt_half          = (real)0.5*dt;
   const real *Passive         = NULL;                // RTVD does not support passive scalars
         bool RuleOut          = false;
   const bool CheckMinPres_Yes = true;

   real   _rho, vx, p, c, Temp, Fluid[5], Fluid_half[5];
   int    ID1, ID2, ID3, Comp[5], delta_k;


// set the order of component for update in different directions
   switch ( XYZ )
   {
      case 0:  Comp[0] = 0;   Comp[1] = 1;   Comp[2] = 2;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 3:  Comp[0] = 0;   Comp[1] = 2;   Comp[2] = 1;   Comp[3] = 3;   Comp[4] = 4;   break;
      case 6:  Comp[0] = 0;   Comp[1] = 3;   Comp[2] = 2;   Comp[3] = 1;   Comp[4] = 4;   break;
   }


// start the TVD scheme
   do
   {

//    determine the array indices for updating in different directions
      switch ( XYZ )
      {
         case 0:  ID1 = to1D1( k, j, i );    break;
         case 3:  ID1 = to1D1( k, i, j );    break;
         case 6:  ID1 = to1D1( i, k, j );    break;
      }


//    load data into per-thread registers
      for (int v=0; v<5; v++)    Fluid[v] = g_Fluid_In[bx][ Comp[v] ][ID1];


//    a. Evaluate the half-step values of fluid variables
//-----------------------------------------------------------------------------

//    (a1). set variables defined in the center of cell
      _rho = (real)1.0 / Fluid[0];
      vx   = _rho * Fluid[1];
      p    = Hydro_Con2Pres( Fluid[0], Fluid[1], Fluid[2], Fluid[3], Fluid[4], Passive,
                             CheckMinPres_Yes, MinPres, NULL_REAL, EoS_DensEint2Pres, NULL, NULL,
                             EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(p) )
         printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 p, __FILE__, __LINE__, __FUNCTION__ );

      if ( Hydro_CheckNegative(Fluid[0]) )
         printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Fluid[0], __FILE__, __LINE__, __FUNCTION__ );
#     endif
      c    = FABS( vx ) + SQRT(  EoS_DensPres2CSqr( Fluid[0], p, Passive, EoS_AuxArray_Flt,
                                                    EoS_AuxArray_Int, EoS_Table )  );

      s_cw[ty][0][i] = Fluid[1];
      s_cw[ty][1][i] = Fluid[1]*vx + p;
      s_cw[ty][2][i] = Fluid[2]*vx;
      s_cw[ty][3][i] = Fluid[3]*vx;
      s_cw[ty][4][i] = ( Fluid[4]+p )*vx;

      s_cu[ty][0][i] = c*Fluid[0];
      s_cu[ty][1][i] = c*Fluid[1];
      s_cu[ty][2][i] = c*Fluid[2];
      s_cu[ty][3][i] = c*Fluid[3];
      s_cu[ty][4][i] = c*Fluid[4];

      __syncthreads();


//    (a2). set flux defined in the right-hand surface of cell by the upwind scheme
      if ( i < FLU_NXT-1 )
      {
         for (int v=0; v<5; v++)
            s_flux[ty][v][i] = (real)0.5*(  ( s_cu[ty][v][i ]+s_cw[ty][v][i ] ) -
                                            ( s_cu[ty][v][ip]-s_cw[ty][v][ip] )  );
      }

      __syncthreads();


//    (a3). evaluate the intermidiate values (u_half)
//    if ( i > 0 )
      if ( i > 0  &&  i < FLU_NXT-1 )
      {
         for (int v=0; v<5; v++)    Fluid_half[v] = Fluid[v] - _dh*dt_half*( s_flux[ty][v][i] - s_flux[ty][v][im] ) ;

//       apply density and internal energy floors
         Fluid_half[0] = FMAX( Fluid_half[0], MinDens );
         Fluid_half[4] = Hydro_CheckMinEintInEngy( Fluid_half[0], Fluid_half[1], Fluid_half[2], Fluid_half[3], Fluid_half[4],
                                                   MinEint, NULL_REAL );
      }



//    Evaluate the full-step values of fluid variables
//-----------------------------------------------------------------------------

//    (b1). reset variables defined in the center of cell at the intermidate state
      if ( i > 0  &&  i < FLU_NXT-1 )
      {
         _rho = (real)1.0 / Fluid_half[0];
         vx   = _rho * Fluid_half[1];
         p    = Hydro_Con2Pres( Fluid_half[0], Fluid_half[1], Fluid_half[2], Fluid_half[3], Fluid_half[4], Passive,
                                CheckMinPres_Yes, MinPres, NULL_REAL, EoS_DensEint2Pres, NULL, NULL,
                                EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(p) )
            printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    p, __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(Fluid_half[0]) )
            printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid_half[0], __FILE__, __LINE__, __FUNCTION__ );
#        endif

         c    = FABS( vx ) + SQRT(  EoS_DensPres2CSqr( Fluid_half[0], p, Passive, EoS_AuxArray_Flt,
                                                       EoS_AuxArray_Int, EoS_Table )  );

         s_cw[ty][0][i] = Fluid_half[1];
         s_cw[ty][1][i] = Fluid_half[1]*vx + p;
         s_cw[ty][2][i] = Fluid_half[2]*vx;
         s_cw[ty][3][i] = Fluid_half[3]*vx;
         s_cw[ty][4][i] = ( Fluid_half[4]+p )*vx;

         s_cu[ty][0][i] = c*Fluid_half[0];
         s_cu[ty][1][i] = c*Fluid_half[1];
         s_cu[ty][2][i] = c*Fluid_half[2];
         s_cu[ty][3][i] = c*Fluid_half[3];
         s_cu[ty][4][i] = c*Fluid_half[4];
      } // if ( i > 0  &&  i < FLU_NXT-1 )


//    (b2). set the right-moving flux defined in the right-hand surface by the TVD scheme
      if ( i > 0  &&  i < FLU_NXT-2 )
      {
         for (int v=0; v<5; v++)    s_RLflux[ty][v][i] = (real)0.5*( s_cu[ty][v][i] + s_cw[ty][v][i] );
      }

      __syncthreads();


      if ( i > 1  &&  i < FLU_NXT-3 )
      {
         for (int v=0; v<5; v++)
         {
            s_flux[ty][v][i] = s_RLflux[ty][v][i];

            Temp =   ( s_RLflux[ty][v][ip]-s_RLflux[ty][v][i ] )
                   * ( s_RLflux[ty][v][i ]-s_RLflux[ty][v][im] );

            if ( Temp > (real)0.0 )
               s_flux[ty][v][i] += Temp / ( s_RLflux[ty][v][ip]-s_RLflux[ty][v][im] );
         }
      }

      __syncthreads();


//    (b3). set the left-moving flux defined in the left-hand surface by the TVD scheme, get the total flux
//    if ( i < FLU_NXT-2 )
      if ( i > 0  &&  i < FLU_NXT-2 )
      {
         for (int v=0; v<5; v++)    s_RLflux[ty][v][i] = (real)0.5*( s_cu[ty][v][ip] - s_cw[ty][v][ip] );
      }

      __syncthreads();


      if ( i > 1  &&  i < FLU_NXT-3 )
      {
         for (int v=0; v<5; v++)
         {
            s_flux[ty][v][i] -= s_RLflux[ty][v][i];

            Temp =   ( s_RLflux[ty][v][im]-s_RLflux[ty][v][i ] )
                   * ( s_RLflux[ty][v][i ]-s_RLflux[ty][v][ip] );

            if ( Temp > (real)0.0 )
               s_flux[ty][v][i] -= Temp / ( s_RLflux[ty][v][im]-s_RLflux[ty][v][ip] );
         }
      }

      __syncthreads();


//    (b4). advance fluid by one full time-step
//    if ( i > 2 )
//    if ( i > 2  &&  i < FLU_NXT-3 )
      if ( i > 2  &&  i < FLU_NXT-3  &&  RuleOut == false )
      {
         for (int v=0; v<5; v++)    Fluid[v] -= _dh*dt*( s_flux[ty][v][i] - s_flux[ty][v][im] ) ;

//       apply density and internal energy floors
         Fluid[0] = FMAX( Fluid[0], MinDens );
         Fluid[4] = Hydro_CheckMinEintInEngy( Fluid[0], Fluid[1], Fluid[2], Fluid[3], Fluid[4],
                                              MinEint, NULL_REAL );


//       check negative density and energy
#        ifdef CHECK_NEGATIVE_IN_FLUID
         if ( Hydro_CheckNegative(Fluid[0]) )
            printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid[0], __FILE__, __LINE__, __FUNCTION__ );

         if ( Hydro_CheckNegative(Fluid[4]) )
            printf( "ERROR : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                    Fluid[4], __FILE__, __LINE__, __FUNCTION__ );
#        endif


//       store the updated data back to the global memory
         if ( FinalOut )
         {
            switch ( XYZ )
            {
               case 0:  ID2 = to1D2( k, j, i );    break;
               case 3:  ID2 = to1D2( k, i, j );    break;
               case 6:  ID2 = to1D2( i, k, j );    break;
            }

            for (int v=0; v<5; v++)    g_Fluid_Out[bx][ Comp[v] ][ID2] = Fluid[v];
         }

         else
            for (int v=0; v<5; v++)    g_Fluid_In [bx][ Comp[v] ][ID1] = Fluid[v];
      }


//    (b5). save the flux required by the flux-correction operation
      if ( StoreFlux )
      if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
      if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
      if ( i == 0 )
      {
         ID3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

         for (int v=0; v<5; v++)
         {
            g_Flux[bx][XYZ+0][v][ID3] = s_flux[ty][ Comp[v] ][          2];
            g_Flux[bx][XYZ+1][v][ID3] = s_flux[ty][ Comp[v] ][FLU_NXT/2-1];
            g_Flux[bx][XYZ+2][v][ID3] = s_flux[ty][ Comp[v] ][FLU_NXT - 4];
         }
      }


//    reset the target array indices
      j += dj;

      if ( j >= j_end )
      {
         delta_k   = ( j - j_end )/size_j + 1;
         k        += delta_k;
         j        -= __umul24( size_j, delta_k );
      }

      Column0 += dj;


//    if the index k exceeds the maximum allowed value --> reset (j,k) to harmless values and wait for other
//    threads (all threads must exist the while loop "at the same time", otherwise __syncthreads will fail !!)
      if ( k >= k_end )
      {
         j       = 0;
         k       = 0;
         RuleOut = true;
      }


      __syncthreads();

   }
   while ( Column0 < NColumn );

} // FUNCTION : CUFLU_Advance



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == RTVD )
