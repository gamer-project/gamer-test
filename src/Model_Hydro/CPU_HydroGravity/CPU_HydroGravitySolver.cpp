#include "CUPOT.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )


// external functions
#ifdef __CUDACC__
 
#include "../GPU_Hydro/CUFLU_Shared_FluUtility.cu"
 
#else
# include "../../../include/Prototype.h"
#endif // #ifdef __CUDACC__ ... else ...


// include c_ExtAcc_AuxArray[]
#ifdef __CUDACC__
#include "CUDA_ConstMemory.h"
#endif




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_HydroGravitySolver
// Description :  1. Advances the momentum and energy density of a group of patches by gravitational acceleration
//                   (including external gravity)
//                2. Note the SRHD still use the Newtonian gravity as a source term in the relativistic Euler equations.
//                   --> Hence SRHD with the Newtonian gravity can only simulate the low-mass relativistic component
//                       embedded in giant gas sphere.
//                       e.g., relativistic AGN jet in gas sphere.
//
// Note        :  1. Currently this function does NOT ensure the consistency between internal energy and
//                   dual-energy variable (e.g., entropy)
//                   --> This consistency breaks only for cells with the dual-energy status labelled
//                       as DE_UPDATED_BY_ETOT_GRA
//                   --> We restore this consistency in Gra_Close()
//                2. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_New   : Array to store the input and output fluid variables
//                g_Pot_Array_New   : Array storing the input potential (at the current step)
//                                    --> _New: to be distinguishable from g_Pot_Array_USG[], which is defined at the previous step
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                g_Pot_Array_USG   : Array storing the input potential          for UNSPLIT_GRAVITY (at the previous step)
//                g_Flu_Array_USG   : Array storing the input density + momentum for UNSPLIT_GRAVITY (at the previous step)
//                g_DE_Array        : Array storing the dual-energy status (for both input and output)
//                g_Emag_Array      : Array storing the cell-centered magnetic energy
//                                    --> Only for checking minimum internal energy in MHD
//                NPatchGroup       : Number of input patch groups (for CPU only)
//                dt                : Time interval to advance solution
//                dh                : Cell size
//                P5_Gradient       : Use 5-points stencil to evaluate the potential gradient
//                UsePot            : Add self-gravity and/or external potential
//                ExtAcc            : Add external acceleration
//                ExtAcc_Func       : Function pointer to the external acceleration routine (for both CPU and GPU)
//                c_ExtAcc_AuxArray : Auxiliary array for adding external acceleration (for CPU only)
//                                    --> When using GPU, this array is stored in the constant memory header
//                                        CUDA_ConstMemory.h and does not need to be passed as a function argument
//                TimeNew           : Physical time at the current  step (for the external gravity solver)
//                TimeOld           : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//                MinEint           : Internal energy floor
//
// Return      :  g_Flu_Array_New, g_DE_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_Emag_Array   [][ CUBE(PS1) ],
   const real dt, const real dh, const bool P5_Gradient,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double TimeNew, const double TimeOld, const real MinEint,
   const EoS_GUESS_t EoS_GuessHTilde_Func,
   const EoS_TEM2H_t EoS_Temp2HTilde_Func,
   const EoS_H2TEM_t EoS_HTilde2Temp_Func )
#else
void CPU_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_Emag_Array   [][ CUBE(PS1) ],
   const int NPatchGroup,
   const real dt, const real dh, const bool P5_Gradient,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const double TimeNew, const double TimeOld, const real MinEint,
   const EoS_GUESS_t EoS_GuessHTilde_Func,
   const EoS_TEM2H_t EoS_Temp2HTilde_Func,
   const EoS_H2TEM_t EoS_HTilde2Temp_Func,
   const double c_EoS_AuxArray_Flt[],
   const int    c_EoS_AuxArray_Int[],
   const real *const c_EoS_Table[EOS_NTABLE_MAX] )
#endif
{

// check
#  ifdef GAMER_DEBUG
   if ( ExtAcc  &&  TimeNew < 0.0 )
      printf( "ERROR : incorrect TimeNew (%14.7e) !!\n", TimeNew );

#  ifdef UNSPLIT_GRAVITY
   if ( g_Flu_Array_USG == NULL )
      printf( "ERROR : g_Flu_Array_USG == NULL !!\n" );

   if ( UsePot  &&  g_Pot_Array_USG == NULL )
      printf( "ERROR : g_Pot_Array_USG == NULL !!\n" );

   if (  ExtAcc  &&  ( TimeOld >= TimeNew || TimeOld < 0.0 )  )
      printf( "ERROR : incorrect time (TimeOld %14.7e, TimeNew = %14.7e) !!\n", TimeOld, TimeNew );
#  endif

#  ifdef DUAL_ENERGY
   if ( g_DE_Array == NULL )
      printf( "ERROR : g_DE_Array == NULL !!\n" );
#  endif

#  ifdef MHD
   if ( g_Emag_Array == NULL )
      printf( "ERROR : g_Emag_Array == NULL !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG

#if ( defined UNSPLIT_GRAVITY && defined SRHD )
# error: SRHD do not support UNSPLIT_GRAVITY !!
#endif

   const real Gra_Const   = ( P5_Gradient ) ? -dt/(12.0*dh) : -dt/(2.0*dh);
   const int  PS1_sqr     = SQR(PS1);
   const int  didx_new[3] = { 1, GRA_NXT,   SQR(GRA_NXT)   };
#  ifdef UNSPLIT_GRAVITY
   const int  didx_old[3] = { 1, USG_NXT_G, SQR(USG_NXT_G) };
#  endif


// load potential from global to shared memory to improve the GPU performance
#  ifdef __CUDACC__
   __shared__ real s_pot_new[ CUBE(GRA_NXT) ];
#  ifdef UNSPLIT_GRAVITY
   __shared__ real s_pot_old[ CUBE(USG_NXT_G) ];
#  endif

   if ( UsePot )
   {
      for (int t=threadIdx.x; t<CUBE(GRA_NXT); t+=GRA_BLOCK_SIZE)
         s_pot_new[t] = g_Pot_Array_New[blockIdx.x][t];

#     ifdef UNSPLIT_GRAVITY
      for (int t=threadIdx.x; t<CUBE(USG_NXT_G); t+=GRA_BLOCK_SIZE)
         s_pot_old[t] = g_Pot_Array_USG[blockIdx.x][t];
#     endif
   }

   __syncthreads();
#  endif // #ifdef __CUDACC__


// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//     to work on different patches
#  ifdef __CUDACC__
   const int P = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatchGroup*8; P++)
#  endif
   {
//    point to the potential array of the target patch
#     ifdef __CUDACC__
      const real *const pot_new = s_pot_new;
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = s_pot_old;
#     endif
#     else // #ifdef __CUDACC__
      const real *const pot_new = g_Pot_Array_New[P];
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = g_Pot_Array_USG[P];
#     endif
#     endif // #ifdef __CUDACC__ ... else ...


//    loop over all cells of the target patch
//    _g0: indices for the arrays without any ghost zone
      CGPU_LOOP( idx_g0, CUBE(PS1) )
      {
#        if ( defined SRHD && defined UNSPLIT_GRAVITY )
         real Etot_in;
#        elif ( !defined SRHD )
//       Enki = non-kinetic energy (i.e. Etot - Ekin)
         real Enki_in, Ekin_out, Etot_in, _rho2;
#        endif

         real acc_new[3]={0.0, 0.0, 0.0}, px_new, py_new, pz_new, rho_new, Etot_out;
#        ifdef UNSPLIT_GRAVITY
         real acc_old[3]={0.0, 0.0, 0.0}, px_old, py_old, pz_old, rho_old, Emag_in=0.0;
#        endif

         const int i_g0    = idx_g0 % PS1;
         const int j_g0    = idx_g0 % PS1_sqr / PS1;
         const int k_g0    = idx_g0 / PS1_sqr;

         const int i_new   = i_g0 + GRA_GHOST_SIZE;
         const int j_new   = j_g0 + GRA_GHOST_SIZE;
         const int k_new   = k_g0 + GRA_GHOST_SIZE;
         const int idx_new = IDX321( i_new, j_new, k_new, GRA_NXT, GRA_NXT );

#        ifdef UNSPLIT_GRAVITY
         const int i_old   = i_g0 + USG_GHOST_SIZE_G;
         const int j_old   = j_g0 + USG_GHOST_SIZE_G;
         const int k_old   = k_g0 + USG_GHOST_SIZE_G;
         const int idx_old = IDX321( i_old, j_old, k_old, USG_NXT_G, USG_NXT_G );
#        endif


//       external acceleration
         if ( ExtAcc )
         {
            double x, y, z;

            x = g_Corner_Array[P][0] + (double)(i_g0*dh);
            y = g_Corner_Array[P][1] + (double)(j_g0*dh);
            z = g_Corner_Array[P][2] + (double)(k_g0*dh);

            ExtAcc_Func( acc_new, x, y, z, TimeNew, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_new[d] *= dt;

#           ifdef UNSPLIT_GRAVITY
            ExtAcc_Func( acc_old, x, y, z, TimeOld, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_old[d] *= dt;
#           endif
         }


//       self-gravity and external potential
         if ( UsePot )
         {
            const int ip1_new = idx_new + didx_new[0];
            const int jp1_new = idx_new + didx_new[1];
            const int kp1_new = idx_new + didx_new[2];
            const int im1_new = idx_new - didx_new[0];
            const int jm1_new = idx_new - didx_new[1];
            const int km1_new = idx_new - didx_new[2];

#           ifdef UNSPLIT_GRAVITY
            const int ip1_old = idx_old + didx_old[0];
            const int jp1_old = idx_old + didx_old[1];
            const int kp1_old = idx_old + didx_old[2];
            const int im1_old = idx_old - didx_old[0];
            const int jm1_old = idx_old - didx_old[1];
            const int km1_old = idx_old - didx_old[2];
#           endif

            if ( P5_Gradient )
            {
               const real Const_8 = (real)8.0;

               const int ip2_new = ip1_new + didx_new[0];
               const int jp2_new = jp1_new + didx_new[1];
               const int kp2_new = kp1_new + didx_new[2];
               const int im2_new = im1_new - didx_new[0];
               const int jm2_new = jm1_new - didx_new[1];
               const int km2_new = km1_new - didx_new[2];

#              ifdef UNSPLIT_GRAVITY
               const int ip2_old = ip1_old + didx_old[0];
               const int jp2_old = jp1_old + didx_old[1];
               const int kp2_old = kp1_old + didx_old[2];
               const int im2_old = im1_old - didx_old[0];
               const int jm2_old = jm1_old - didx_old[1];
               const int km2_old = km1_old - didx_old[2];
#              endif

               acc_new[0] += Gra_Const*( - pot_new[ip2_new] + Const_8*pot_new[ip1_new] - Const_8*pot_new[im1_new] + pot_new[im2_new] );
               acc_new[1] += Gra_Const*( - pot_new[jp2_new] + Const_8*pot_new[jp1_new] - Const_8*pot_new[jm1_new] + pot_new[jm2_new] );
               acc_new[2] += Gra_Const*( - pot_new[kp2_new] + Const_8*pot_new[kp1_new] - Const_8*pot_new[km1_new] + pot_new[km2_new] );

#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( - pot_old[ip2_old] + Const_8*pot_old[ip1_old] - Const_8*pot_old[im1_old] + pot_old[im2_old] );
               acc_old[1] += Gra_Const*( - pot_old[jp2_old] + Const_8*pot_old[jp1_old] - Const_8*pot_old[jm1_old] + pot_old[jm2_old] );
               acc_old[2] += Gra_Const*( - pot_old[kp2_old] + Const_8*pot_old[kp1_old] - Const_8*pot_old[km1_old] + pot_old[km2_old] );
#              endif
            } // if ( P5_Gradient )

            else
            {
               acc_new[0] += Gra_Const*( pot_new[ip1_new] - pot_new[im1_new] );
               acc_new[1] += Gra_Const*( pot_new[jp1_new] - pot_new[jm1_new] );
               acc_new[2] += Gra_Const*( pot_new[kp1_new] - pot_new[km1_new] );

#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( pot_old[ip1_old] - pot_old[im1_old] );
               acc_old[1] += Gra_Const*( pot_old[jp1_old] - pot_old[jm1_old] );
               acc_old[2] += Gra_Const*( pot_old[kp1_old] - pot_old[km1_old] );
#              endif
            } // if ( P5_Gradient ) ... else ...
         } // if ( UsePot )

#        ifdef SRHD
         real Cons_new[NCOMP_FLUID], Prim_new[NCOMP_FLUID], LorentzFactor_new;

         Cons_new[DENS] = g_Flu_Array_New[P][DENS][idx_g0];
         Cons_new[MOMX] = g_Flu_Array_New[P][MOMX][idx_g0];
         Cons_new[MOMY] = g_Flu_Array_New[P][MOMY][idx_g0];
         Cons_new[MOMZ] = g_Flu_Array_New[P][MOMZ][idx_g0];
         Cons_new[ENGY] = g_Flu_Array_New[P][ENGY][idx_g0];

         Hydro_Con2Pri( Cons_new, Prim_new, NULL_REAL, NULL_BOOL, NULL_INT, NULL,
                        NULL_BOOL, NULL_REAL, NULL, NULL, EoS_GuessHTilde_Func,
                        EoS_HTilde2Temp_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int,
                        c_EoS_Table, NULL, &LorentzFactor_new );

         rho_new = Prim_new[0];
         px_new  = Cons_new[MOMX];
         py_new  = Cons_new[MOMY];
         pz_new  = Cons_new[MOMZ];

#        ifdef UNSPLIT_GRAVITY
         real Vx_new, Vy_new, Vz_new;
         Vx_new  = Prim_new[1] / LorentzFactor_new;
         Vy_new  = Prim_new[2] / LorentzFactor_new;
         Vz_new  = Prim_new[3] / LorentzFactor_new;
#        endif
#        endif

#        if ( ( defined SRHD && defined UNSPLIT_GRAVITY ) || !defined SRHD )
         Etot_in = g_Flu_Array_New[P][ENGY][idx_g0];
#        endif

//       advance fluid
#        ifdef UNSPLIT_GRAVITY


#        ifdef SRHD
         real Cons_old[NCOMP_FLUID], Prim_old[NCOMP_FLUID], LorentzFactor_old;
         real Vx_old, Vy_old, Vz_old;

         Cons_old[DENS] = g_Flu_Array_USG[P][DENS][idx_g0];
         Cons_old[MOMX] = g_Flu_Array_USG[P][MOMX][idx_g0];
         Cons_old[MOMY] = g_Flu_Array_USG[P][MOMY][idx_g0];
         Cons_old[MOMZ] = g_Flu_Array_USG[P][MOMZ][idx_g0];
         Cons_old[ENGY] = g_Flu_Array_USG[P][ENGY][idx_g0];

         Hydro_Con2Pri( Cons_old, Prim_old, NULL_REAL, NULL_BOOL, NULL_INT, NULL,
                        NULL_BOOL, NULL_REAL, NULL, NULL, EoS_GuessHTilde_Func,
                        EoS_HTilde2Temp_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int,
                        c_EoS_Table, NULL, &LorentzFactor_old );
 
         rho_old = Prim_old[0];
         Vx_old  = Prim_old[1] / LorentzFactor_old;
         Vy_old  = Prim_old[2] / LorentzFactor_old;
         Vz_old  = Prim_old[3] / LorentzFactor_old;
         px_old  = Cons_old[MOMX];
         py_old  = Cons_old[MOMY];
         pz_old  = Cons_old[MOMZ];
 
#        else
         rho_new = g_Flu_Array_New[P][DENS][idx_g0];
         rho_old = g_Flu_Array_USG[P][DENS][idx_g0];
         px_new  = g_Flu_Array_New[P][MOMX][idx_g0];
         px_old  = g_Flu_Array_USG[P][MOMX][idx_g0];
         py_new  = g_Flu_Array_New[P][MOMY][idx_g0];
         py_old  = g_Flu_Array_USG[P][MOMY][idx_g0];
         pz_new  = g_Flu_Array_New[P][MOMZ][idx_g0];
         pz_old  = g_Flu_Array_USG[P][MOMZ][idx_g0];

//       backup the original non-kinetic energy so that we can restore it later if necessary
         _rho2   = (real)0.5/rho_new;
         Enki_in = Etot_in - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
#        ifdef MHD
         Emag_in = g_Emag_Array[P][idx_g0];
#        endif // # ifdef MHD
#        endif // # ifdef SRHD

//       update the momentum density
         px_new += (real)0.5*( rho_old*acc_old[0] + rho_new*acc_new[0] );
         py_new += (real)0.5*( rho_old*acc_old[1] + rho_new*acc_new[1] );
         pz_new += (real)0.5*( rho_old*acc_old[2] + rho_new*acc_new[2] );

         g_Flu_Array_New[P][MOMX][idx_g0] = px_new;
         g_Flu_Array_New[P][MOMY][idx_g0] = py_new;
         g_Flu_Array_New[P][MOMZ][idx_g0] = pz_new;

#        ifndef SRHD
//       record the updated kinematic energy density
         Ekin_out = _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
#        endif

//       update the total energy density
#        ifdef DUAL_ENERGY

//       for the unsplitting method with the dual-energy formalism, we correct the **total energy density**
//       only if the dual-energy status != DE_UPDATED_BY_DUAL
//       --> for (a) DE_UPDATED_BY_DUAL     --> Eint has been updated by the dual-energy variable
//               (b) DE_UPDATED_BY_MIN_PRES --> Eint has been set to the minimum threshold
//       --> currently for (b) we still update the total energy density

         if ( g_DE_Array[P][idx_g0] == DE_UPDATED_BY_DUAL )
         {
//          fix the internal energy and the dual-energy variable
            Etot_out = Enki_in + Ekin_out;
         }

         else
         {
//          update the total energy, where internal energy and dual-energy variable may change as well
            Etot_out = Etot_in + (real)0.5*( px_old*acc_old[0] + py_old*acc_old[1] + pz_old*acc_old[2] +
                                             px_new*acc_new[0] + py_new*acc_new[1] + pz_new*acc_new[2] );

//          check the minimum internal energy
//###NOTE: assuming Etot = Eint + Ekin + Emag
//          (a) if the updated internal energy is greater than the threshold, set the dual-energy status == DE_UPDATED_BY_ETOT_GRA
            if ( Etot_out - Ekin_out - Emag_in >= MinEint )
               g_DE_Array[P][idx_g0] = DE_UPDATED_BY_ETOT_GRA;

//          (b) otherwise restore the original internal energy and keep the original dual-energy status
            else
               Etot_out = Enki_in + Ekin_out;
         }

#        else // # ifdef DUAL_ENERGY

#        ifdef SRHD
         px_new = rho_new*Vx_new;
         py_new = rho_new*Vy_new; 
         pz_new = rho_new*Vz_new; 
         px_old = rho_old*Vx_old; 
         py_old = rho_old*Vy_old; 
         pz_old = rho_old*Vz_old; 
#        endif

//       for the unsplitting method without the dual-energy formalism, we always correct the total energy density
//       instead of the kinematic energy density
//       --> internal energy may change
//       --> we must check the minimum internal energy after this update
         Etot_out = Etot_in + (real)0.5*( px_old*acc_old[0] + py_old*acc_old[1] + pz_old*acc_old[2] +
                                          px_new*acc_new[0] + py_new*acc_new[1] + pz_new*acc_new[2] );

#        ifndef SRHD
//       check the minimum internal energy
//       --> restore the original internal energy if the updated value becomes smaller than the threshold
         if ( Etot_out - Ekin_out - Emag_in < MinEint )
            Etot_out = Enki_in + Ekin_out;
#        endif

#        endif // #ifdef DUAL_ENERGY ... else ...


#        else  // #ifdef UNSPLIT_GRAVITY

#        ifdef SRHD
//       backup the dimensionless temperature (k_{B}T/mc^2) so that we can restore it later
         real Temperature = Prim_new[4]/Prim_new[0];
#        else
         rho_new = g_Flu_Array_New[P][DENS][idx_g0];
         px_new  = g_Flu_Array_New[P][MOMX][idx_g0];
         py_new  = g_Flu_Array_New[P][MOMY][idx_g0];
         pz_new  = g_Flu_Array_New[P][MOMZ][idx_g0];

//       backup the original internal energy so that we can restore it later
         _rho2   = (real)0.5/rho_new;
         Enki_in = Etot_in - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
#        endif

//       update the momentum density
         px_new += rho_new*acc_new[0];
         py_new += rho_new*acc_new[1];
         pz_new += rho_new*acc_new[2];

         g_Flu_Array_New[P][MOMX][idx_g0] = px_new;
         g_Flu_Array_New[P][MOMY][idx_g0] = py_new;
         g_Flu_Array_New[P][MOMZ][idx_g0] = pz_new;

//       for the splitting method, we ensure that the internal energy is unchanged
#        ifdef SRHD
         real Msqr = SQR(px_new) + SQR(py_new) + SQR(pz_new);
         real Dsqr = SQR(Cons_new[DENS]);
         real HTilde = EoS_Temp2HTilde_Func( Temperature, NULL, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table );
         real h = HTilde + (real)1.0;
         real factor = SQRT( Dsqr*h*h + Msqr );

         Cons_new[ENGY]  = Dsqr * SQR(HTilde) + (real)2.0*HTilde*Dsqr + Msqr;
         Cons_new[ENGY] /= ( factor + Cons_new[DENS] );
         Cons_new[ENGY] -= Dsqr*h*Temperature / factor;
         Etot_out = Cons_new[ENGY];
#        else
         Ekin_out = _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
         Etot_out = Enki_in + Ekin_out;
#        endif


#        endif // #ifdef UNSPLIT_GRAVITY ... else ...



//       store the updated total energy density to the output array
         g_Flu_Array_New[P][ENGY][idx_g0] = Etot_out;


#        ifdef SRHD
#        ifdef CHECK_FAILED_CELL_IN_FLUID
//       check unphysical result
         real Cons[NCOMP_FLUID];

         Cons[DENS] = g_Flu_Array_New[P][DENS][idx_g0];
         Cons[MOMX] = g_Flu_Array_New[P][MOMX][idx_g0];
         Cons[MOMY] = g_Flu_Array_New[P][MOMY][idx_g0];
         Cons[MOMZ] = g_Flu_Array_New[P][MOMZ][idx_g0];
         Cons[ENGY] = g_Flu_Array_New[P][ENGY][idx_g0];

         SRHydro_CheckUnphysical( Cons, NULL, __FUNCTION__, __LINE__, true );
#        endif
#        endif

      } // CGPU_LOOP( idx_g0, CUBE(PS1) )
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_HydroGravitySolver



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
