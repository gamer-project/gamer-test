#include "GAMER.h"
#if ( MODEL == HYDRO )


// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Function_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] );

// this function pointer must be set by a test problem initializer
void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                const int lv, double AuxArray[] ) = NULL;

extern bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] );

#ifdef MHD
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_Function_BField_User_Template( real magnetic[], const double x, const double y, const double z, const double Time,
                                                const int lv, double AuxArray[] );

// this function pointer must be set by a test problem initializer
void (*Init_Function_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] ) = NULL;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User_Template
// Description :  Function template to initialize the fluid field
//
// Note        :  1. Invoked by Hydro_Init_ByFunction_AssignData() using the function pointer
//                   "Init_Function_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in Hydro_Init_ByFunction_AssignData()
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User_Template( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] )
{

   const real Dens0 = 1.0;
   const real Vx0   = 1.25e-1;
   const real Vy0   = 2.30e-1;
   const real Vz0   = 3.70e-1;
   const real Pres0 = 1.0;
   const real Emag0 = 0.0;    // must be zero here even for MHD

   real Dens, Vx, Vy, Vz, Pres;
   real MomX, MomY, MomZ, Eint, Etot;
#  if ( NCOMP_PASSIVE > 0 )
   real Passive[NCOMP_PASSIVE];
#  else
   real *Passive = NULL;
#  endif

// set the primitive variables
   Dens = Dens0 + 0.2*exp(  -(  SQR(1.1*x-0.5*amr->BoxSize[0])
                               +SQR(2.2*y-0.5*amr->BoxSize[1])
                               +SQR(3.3*z-0.5*amr->BoxSize[2]) ) / SQR( 1.8*amr->BoxSize[2] )  );
   Vx   = Vx0*sin( 2.0*M_PI/amr->BoxSize[0] );
   Vy   = Vy0*cos( 2.0*M_PI/amr->BoxSize[1] );
   Vz   = Vz0*sin( 2.0*M_PI/amr->BoxSize[2] );
   Pres = Pres0*(  2.0 + sin( 2.0*M_PI*(4.5*x+5.5*y*6.5*z)/amr->BoxSize[2] )  );

// set passive scalars
#  if ( NCOMP_PASSIVE > 0 )
// Passive[X] = ...;
#  endif

// convert primitive variables to conservative variables
   MomX = Dens*Vx;
   MomY = Dens*Vy;
   MomZ = Dens*Vz;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, Emag0 );

// store the results
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
#  if ( NCOMP_PASSIVE > 0 )
// fluid[XXXX] = ...;
#  endif

} // FUNCTION : Init_Function_User_Template



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_BField_User_Template
// Description :  Function template to initialize the magnetic field
//
// Note        :  1. Invoked by Hydro_Init_ByFunction_AssignData() using the function pointer
//                   "Init_Function_BField_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (uless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void Init_Function_BField_User_Template( real magnetic[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] )
{

   magnetic[MAGX] = 1.0;
   magnetic[MAGY] = 2.0;
   magnetic[MAGZ] = 3.0;

} // FUNCTION : Init_Function_BField_User_Template
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  HydroInit_ByFunction_AssignData
// Description :  Construct the initial condition in HYDRO
//
// Note        :  1. Work for the option "OPT__INIT == INIT_BY_FUNCTION"
//                2. Function pointers "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr/Init_Function_BField_User_Ptr"
//                   must be set by a test problem initializer
//                3. One can disable OpenMP in this routine by setting OPT__INIT_GRID_WITH_OMP = 0
//                   --> Useful if "Init_Function_User_Ptr/Flu_ResetByUser_Func_Ptr/Init_Function_BField_User_Ptr"
//                       do not support OpenMP
//                       (e.g., they may not be thread-safe or may involve a random number generator for which
//                       all threads would share the same random seed when adopting OpenMP)
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_ByFunction_AssignData( const int lv )
{

// check
   if ( Init_Function_User_Ptr == NULL )  Aux_Error( ERROR_INFO, "Init_Function_User_Ptr == NULL !!\n" );

#  ifdef MHD
   if ( Init_Function_BField_User_Ptr == NULL  &&  !OPT__INIT_BFIELD_BYFILE )
      Aux_Error( ERROR_INFO, "Init_Function_BField_User_Ptr == NULL !!\n" );
#  endif

   if ( OPT__RESET_FLUID  &&  Flu_ResetByUser_Func_Ptr == NULL )
      Aux_Error( ERROR_INFO, "Flu_ResetByUser_Func_Ptr == NULL for OPT__RESET_FLUID !!\n" );


// set the number of OpenMP threads
#  ifdef OPENMP
   const int OMP_NThread = ( OPT__INIT_GRID_WITH_OMP ) ? OMP_NTHREAD : 1;
#  endif

   const int    NSub     = ( INIT_SUBSAMPLING_NCELL <= 0 ) ? 1 : INIT_SUBSAMPLING_NCELL;
   const double dh       = amr->dh[lv];
   const double dh_sub   = dh / NSub;
   const double _NSub3   = 1.0/CUBE(NSub);
#  ifdef MHD
   const double _NSub2   = 1.0/SQR(NSub);
#  endif


#  ifdef MHD
   if ( OPT__INIT_BFIELD_BYFILE )  MHD_Init_BField_ByFile(lv);
#  endif

#  pragma omp parallel for schedule( runtime ) num_threads( OMP_NThread )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    1. set the magnetic field
#     ifdef MHD

      if ( !OPT__INIT_BFIELD_BYFILE ) {

         real magnetic_1v, magnetic_sub[NCOMP_MAG];

   //    loop over B_X/Y/Z to set one component at a time
   //    --> because different components are defined at different cell faces
         for (int v=0; v<NCOMP_MAG; v++)
         {
            int    ijk_end[3], sub_end[3], idx=0;
            double dxyz0[3];

            for (int d=0; d<3; d++)
            {
               ijk_end[d] = ( d == v ) ? PS1+1 : PS1;
               sub_end[d] = ( d == v ) ? 1     : NSub;
               dxyz0  [d] = ( d == v ) ? 0.0   : 0.5*dh_sub;
            }

            for (int k=0; k<ijk_end[2]; k++)    {  const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + dxyz0[2];
            for (int j=0; j<ijk_end[1]; j++)    {  const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + dxyz0[1];
            for (int i=0; i<ijk_end[0]; i++)    {  const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + dxyz0[0];

               magnetic_1v = (real)0.0;

               for (int kk=0; kk<sub_end[2]; kk++)    {  const double z = z0 + kk*dh_sub;
               for (int jj=0; jj<sub_end[1]; jj++)    {  const double y = y0 + jj*dh_sub;
               for (int ii=0; ii<sub_end[0]; ii++)    {  const double x = x0 + ii*dh_sub;

                  Init_Function_BField_User_Ptr( magnetic_sub, x, y, z, Time[lv], lv, NULL );

                  magnetic_1v += magnetic_sub[v];

               }}}

               amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[v][ idx ++ ] = magnetic_1v*_NSub2;
            }}} // i,j,k
         } // for (int v=0; v<NCOMP_MAG; v++)
      } // if ( !OPT__INIT_BFIELD_BY_FILE )

#     endif // #ifdef MHD



//    2. set the fluid field
      real fluid[NCOMP_TOTAL], fluid_sub[NCOMP_TOTAL];

      for (int k=0; k<PS1; k++)  {  const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh_sub;
      for (int j=0; j<PS1; j++)  {  const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh_sub;
      for (int i=0; i<PS1; i++)  {  const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh_sub;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = (real)0.0;

         for (int kk=0; kk<NSub; kk++)    {  const double z = z0 + kk*dh_sub;
         for (int jj=0; jj<NSub; jj++)    {  const double y = y0 + jj*dh_sub;
         for (int ii=0; ii<NSub; ii++)    {  const double x = x0 + ii*dh_sub;

            Init_Function_User_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

//          modify the initial condition if required
            if ( OPT__RESET_FLUID )
               Flu_ResetByUser_Func_Ptr( fluid_sub, x, y, z, Time[lv], lv, NULL );

            for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] += fluid_sub[v];

         }}}

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] *= _NSub3;


//       add the magnetic energy
#        ifdef MHD
         const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
         fluid[ENGY] += Emag;
#        else
         const real Emag = NULL_REAL;
#        endif


//       apply density and internal energy floors
#        if ( !defined SRHD )
         fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
         fluid[ENGY] = Hydro_CheckMinEintInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                 MIN_EINT, Emag );
#        endif

//       calculate the dual-energy variable (entropy or internal energy)
#        if   ( DUAL_ENERGY == DE_ENPY )
         fluid[ENPY] = Hydro_Con2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Emag,
                                          EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#        elif ( DUAL_ENERGY == DE_EINT )
#        error : DE_EINT is NOT supported yet !!
#        endif

//       floor and normalize passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

         if ( OPT__NORMALIZE_PASSIVE )
            Hydro_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#        endif

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // FUNCTION : Hydro_Init_ByFunction_AssignData



#endif // #if ( MODEL == HYDRO )
