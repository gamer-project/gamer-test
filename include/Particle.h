#ifndef __PARTICLE_H__
#define __PARTICLE_H__



#ifndef PARTICLE
#  error : ERROR : PARTICLE is not defined !!
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_PARTICLE
#endif

// Factors used in "Add(One)Particle" and "Remove(One)Particle" to resize the allocated arrays
// PARLIST_GROWTH_FACTOR must >= 1.0; PARLIST_REDUCE_FACTOR must <= 1.0
#  define PARLIST_GROWTH_FACTOR     1.1
#  define PARLIST_REDUCE_FACTOR     0.8

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  Particle_t
// Description :  Data structure of particles
//
// Data Member :  ParListSize             : Size of the particle data arrays (must be >= NPar_AcPlusInac)
//                InactiveParListSize     : Size of the inactive particle list (InactiveParList)
//                NPar_Active_AllRank     : Total number of active particles summed up over all MPI ranks
//                NPar_AcPlusInac         : Total number of particles (active + inactive) in this MPI rank
//                NPar_Active             : Total number of active particles in this MPI rank
//                                          --> Inactive particles: particles removed from simulations because of
//                                           (1) lying outside the active region (mass == PAR_INACTIVE_OUTSIDE)
//                                           (2) being sent to other MPI ranks   (mass == PAR_INACTIVE_MPI)
//                NPar_Inactive           : Total number of inactive particles in this MPI rank
//                NPar_Lv                 : Total number of active particles at each level in this MPI rank
//                Init                    : Initialization methods (1/2/3 --> call function/restart/load from file)
//                ParICFormat             : Data format of the particle initialization file (1=[att][id], 2=[id][att])
//                ParICMass               : Assign this mass to all particles for Init=3
//                Interp                  : Mass/acceleration interpolation scheme (NGP,CIC,TSC)
//                InterpTracer            : Velocity interpolation scheme (NGP,CIC,TSC)
//                Integ                   : Integration scheme (PAR_INTEG_EULER, PAR_INTEG_KDK)
//                IntegTracer             : Integration scheme for tracers (PAR_INTEG_EULER, PAR_INTEG_RK2)
//                ImproveAcc              : Improve force accuracy around the patch boundaries
//                                          (by using potential in the patch ghost zone instead of nearby patch
//                                          or interpolation)
//                PredictPos              : Predict particle position during mass assignment
//                RemoveCell              : remove particles RemoveCell-base-level-cells away from the boundary
//                                          (for non-periodic BC only)
//                GhostSize               : Number of ghost zones required for interpolation scheme
//                Attribute               : Pointer arrays to different particle attributes (Mass, Pos, Vel, ...)
//                InactiveParList         : List of inactive particle IDs
//                R2B_Real_NPatchTotal    : see R2B_Buff_NPatchTotal
//                R2B_Real_NPatchEachRank : see R2B_Buff_NPatchEachRank
//                R2B_Real_PIDList        : see R2B_Buff_PIDList
//                R2B_Buff_NPatchTotal    : Number of buffer patches to receive particles from the corresponding real patches
//                                          --> R2B : Real to Buffer
//                                          --> Mainly for the Poisson solver
//                                          --> For LOAD_BALANCE only
//                R2B_Buff_NPatchEachRank :  Number of buffer patches to receive data from different MPI ranks
//                R2B_Buff_PIDList        : Buffer patch IDs list to receive particles from other ranks
//                                          --> Mainly for the Poisson solver
//                                          --> For LOAD_BALANCE only
//                                          --> It has the dimension [NLEVEL][2]
//                                              [lv][0/1] is for receiving particles at lv/lv-1 around real patches at lv
//                                              --> Mainly for the Poisson solver at lv (i.e., for calculating the
//                                                  total density field at lv)
//                                              --> More specific,
//                                                  [lv][0] is for receiving the particles of sibling-buffer patches at lv
//                                                  adjacent to real patches at lv
//                                                  [lv][1] is for receiving the particles of father-sibling-buffer patches at lv-1
//                                                  adjacent to real patches at lv
//                                              --> Both are constructed by calling Par_LB_RecordExchangeParticlePatchID( lv )
//                B2R_Real_NPatchTotal    : Similar to R2B_Real_NPatchTotal,    but for sending particles from buffer to real patches
//                B2R_Real_NPatchEachRank : Similar to R2B_Real_NPatchEachRank, ...
//                B2R_Real_PIDList        : Similar to R2B_Real_PIDList,        ...
//                B2R_Buff_NPatchTotal    : Similar to R2B_Buff_NPatchTotal,    ...
//                B2R_Buff_NPatchEachRank : Similar to R2B_Buff_NPatchEachRank, ...
//                B2R_Buff_PIDList        : Similar to R2B_Buff_PIDList,        ...
//                F2S_Send_NPatchTotal    : Similar to R2B_Real_NPatchTotal,    but for sending particles from father to son patches
//                F2S_Send_NPatchEachRank : Similar to R2B_Real_NPatchEachRank, ...
//                F2S_Send_PIDList        : Similar to R2B_Real_PIDList,        ...
//                F2S_Recv_NPatchTotal    : Similar to R2B_Buff_NPatchTotal,    ...
//                F2S_Recv_NPatchEachRank : Similar to R2B_Buff_NPatchEachRank, ...
//                F2S_Recv_PIDList        : Similar to R2B_Buff_PIDList,        ...
//                Mass                    : Particle mass
//                                          Mass < 0.0 --> this particle has been removed from simulations
//                                                     --> PAR_INACTIVE_OUTSIDE: fly outside the simulation box
//                                                         PAR_INACTIVE_MPI:     sent to other MPI ranks
//                Pos                     : Particle position
//                Vel                     : Particle velocity
//                Time                    : Particle physical time
//                Acc                     : Particle acceleration (only when STORE_PAR_ACC is on)
//
// Method      :  Particle_t        : Constructor
//               ~Particle_t        : Destructor
//                InitRepo          : Initialize particle repository
//                AddOneParticle    : Add one new particle into the particle list
//                RemoveOneParticle : Remove one particle from the particle list
//-------------------------------------------------------------------------------------------------------
struct Particle_t
{

// data members
// ===================================================================================
   long          ParListSize;
   long          InactiveParListSize;
   long          NPar_Active_AllRank;
   long          NPar_AcPlusInac;
   long          NPar_Active;
   long          NPar_Inactive;
   long          NPar_Lv[NLEVEL];
   ParInit_t     Init;
   ParICFormat_t ParICFormat;
   double        ParICMass;
   ParInterp_t   Interp;
   ParInteg_t    Integ;
   TracerInteg_t IntegTracer;
   bool          ImproveAcc;
   bool          PredictPos;
   double        RemoveCell;
   int           GhostSize;
   real         *Attribute[PAR_NATT_TOTAL];
   long         *InactiveParList;

#  ifdef LOAD_BALANCE
   int           R2B_Real_NPatchTotal   [NLEVEL][2];
   int          *R2B_Real_NPatchEachRank[NLEVEL][2];
   int          *R2B_Real_PIDList       [NLEVEL][2];
   int           R2B_Buff_NPatchTotal   [NLEVEL][2];
   int          *R2B_Buff_NPatchEachRank[NLEVEL][2];
   int          *R2B_Buff_PIDList       [NLEVEL][2];

   int           B2R_Real_NPatchTotal   [NLEVEL][2];
   int          *B2R_Real_NPatchEachRank[NLEVEL][2];
   int          *B2R_Real_PIDList       [NLEVEL][2];
   int           B2R_Buff_NPatchTotal   [NLEVEL][2];
   int          *B2R_Buff_NPatchEachRank[NLEVEL][2];
   int          *B2R_Buff_PIDList       [NLEVEL][2];

   int           F2S_Send_NPatchTotal   [NLEVEL];
   int          *F2S_Send_NPatchEachRank[NLEVEL];
   int          *F2S_Send_PIDList       [NLEVEL];
   int           F2S_Recv_NPatchTotal   [NLEVEL];
   int          *F2S_Recv_NPatchEachRank[NLEVEL];
   int          *F2S_Recv_PIDList       [NLEVEL];
#  endif // #ifdef LOAD_BALANCE

   real         *Mass;
   real         *PosX;
   real         *PosY;
   real         *PosZ;
   real         *VelX;
   real         *VelY;
   real         *VelZ;
   real         *Time;
   real         *Type;
#  ifdef STORE_PAR_ACC
   real         *AccX;
   real         *AccY;
   real         *AccZ;
#  endif


   //===================================================================================
   // Constructor :  Particle_t
   // Description :  Constructor of the structure "Particle_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  None
   //===================================================================================
   Particle_t()
   {

      NPar_Active_AllRank = -1;
      NPar_AcPlusInac     = -1;
      Init                = PAR_INIT_NONE;
      ParICFormat         = PAR_IC_FORMAT_NONE;
      ParICMass           = -1.0;
      Interp              = PAR_INTERP_NONE;
      Integ               = PAR_INTEG_NONE;
      IntegTracer         = TRACER_INTEG_NONE;
      ImproveAcc          = true;
      PredictPos          = true;
      RemoveCell          = -999.9;
      GhostSize           = -1;

      for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv[lv] = 0;

      for (int v=0; v<PAR_NATT_TOTAL; v++)   Attribute[v] = NULL;

      InactiveParList = NULL;

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
         R2B_Real_NPatchTotal   [lv][t] = 0;
         R2B_Real_NPatchEachRank[lv][t] = NULL;
         R2B_Real_PIDList       [lv][t] = NULL;
         R2B_Buff_NPatchTotal   [lv][t] = 0;
         R2B_Buff_NPatchEachRank[lv][t] = NULL;
         R2B_Buff_PIDList       [lv][t] = NULL;

         B2R_Real_NPatchTotal   [lv][t] = 0;
         B2R_Real_NPatchEachRank[lv][t] = NULL;
         B2R_Real_PIDList       [lv][t] = NULL;
         B2R_Buff_NPatchTotal   [lv][t] = 0;
         B2R_Buff_NPatchEachRank[lv][t] = NULL;
         B2R_Buff_PIDList       [lv][t] = NULL;
         }

         F2S_Send_NPatchTotal   [lv]    = 0;
         F2S_Send_NPatchEachRank[lv]    = NULL;
         F2S_Send_PIDList       [lv]    = NULL;
         F2S_Recv_NPatchTotal   [lv]    = 0;
         F2S_Recv_NPatchEachRank[lv]    = NULL;
         F2S_Recv_PIDList       [lv]    = NULL;
      } // for (int lv=0; lv<NLEVEL; lv++)
#     endif // #ifdef LOAD_BALANCE

      Mass = NULL;
      PosX = NULL;
      PosY = NULL;
      PosZ = NULL;
      VelX = NULL;
      VelY = NULL;
      VelZ = NULL;
      Time = NULL;
      Type = NULL;
#     ifdef STORE_PAR_ACC
      AccX = NULL;
      AccY = NULL;
      AccZ = NULL;
#     endif

   } // METHOD : Particle_t



   //===================================================================================
   // Destructor  :  ~Particle_t
   // Description :  Destructor of the structure "Particle_t"
   //
   // Note        :  Free memory
   //===================================================================================
   ~Particle_t()
   {

      for (int v=0; v<PAR_NATT_TOTAL; v++)
         if ( Attribute[v] != NULL )   free( Attribute[v] );

      if ( InactiveParList != NULL )   free( InactiveParList );

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
         if ( R2B_Real_NPatchEachRank[lv][t] != NULL )   delete [] R2B_Real_NPatchEachRank[lv][t];
         if ( R2B_Buff_NPatchEachRank[lv][t] != NULL )   delete [] R2B_Buff_NPatchEachRank[lv][t];
         if ( R2B_Real_PIDList       [lv][t] != NULL )   delete [] R2B_Real_PIDList       [lv][t];
         if ( R2B_Buff_PIDList       [lv][t] != NULL )   free(     R2B_Buff_PIDList       [lv][t] );

         if ( B2R_Real_NPatchEachRank[lv][t] != NULL )   delete [] B2R_Real_NPatchEachRank[lv][t];
         if ( B2R_Buff_NPatchEachRank[lv][t] != NULL )   delete [] B2R_Buff_NPatchEachRank[lv][t];
         if ( B2R_Real_PIDList       [lv][t] != NULL )   delete [] B2R_Real_PIDList       [lv][t];
         if ( B2R_Buff_PIDList       [lv][t] != NULL )   free(     B2R_Buff_PIDList       [lv][t] );
         } // for (int t=0; t<2; t++)

         if ( F2S_Send_NPatchEachRank[lv] != NULL )      delete [] F2S_Send_NPatchEachRank[lv];
         if ( F2S_Recv_NPatchEachRank[lv] != NULL )      delete [] F2S_Recv_NPatchEachRank[lv];
         if ( F2S_Send_PIDList       [lv] != NULL )      free(     F2S_Send_PIDList       [lv] );
         if ( F2S_Recv_PIDList       [lv] != NULL )      delete [] F2S_Recv_PIDList       [lv];
      } // for (int lv=0; lv<NLEVE; lv++)
#     endif // LOAD_BALANCE

   } // METHOD : ~Particle_t



   //===================================================================================
   // Constructor :  InitRepo
   // Description :  Initialize particle repository
   //                --> All variables related to the number of particles
   //
   // Note        :  1. Initialize both "NPar_AcPlusInac, NPar_Active and ParListSize" as NPar_Input
   //                   --> Assuming no inactive particles (i.e., NPar_Inactive = 0)
   //                2. For LOAD_BALANCE, some lists recording the information for exchanging
   //                   particles between different ranks are also allocated here
   //
   // Parameter   :  NPar_Input : Total number of active particles
   //                NRank      : Total number of MPI ranks
   //===================================================================================
   void InitRepo( const long NPar_Input, const int NRank )
   {

//    check
      if ( NPar_Input < 0 )   Aux_Error( ERROR_INFO, "NPar_Input (%ld) < 0 !!\n", NPar_Input );

//    initialize variables related to the number of particles
      NPar_AcPlusInac     = NPar_Input;
      NPar_Active         = NPar_Input;                  // assuming all particles are active initially
      NPar_Inactive       = 0;
      ParListSize         = NPar_Input;                  // set ParListSize = NPar_AcPlusInac at the beginning
      InactiveParListSize = MAX( 1, ParListSize/100 );   // set arbitrarily (but must > 0)

//    allocate arrays (use malloc so that realloc can be used later to resize the array)
//    --> free memory first since other functions (e.g., LB_Init_LoadBalance()) will call InitRepo() again
      for (int v=0; v<PAR_NATT_TOTAL; v++)
      {
         if ( Attribute[v] != NULL )   free( Attribute[v] );
         Attribute[v] = (real*)malloc( ParListSize*sizeof(real) );
      }

      if ( InactiveParList != NULL )   free( InactiveParList );
      InactiveParList = (long*)malloc( InactiveParListSize*sizeof(long) );

#     ifdef LOAD_BALANCE
      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int t=0; t<2; t++) {
            if ( R2B_Real_NPatchEachRank[lv][t] == NULL )   R2B_Real_NPatchEachRank[lv][t] = new int [NRank];
            if ( R2B_Buff_NPatchEachRank[lv][t] == NULL )   R2B_Buff_NPatchEachRank[lv][t] = new int [NRank];

            for (int r=0; r<NRank; r++)
            {
               R2B_Real_NPatchEachRank[lv][t][r] = 0;
               R2B_Buff_NPatchEachRank[lv][t][r] = 0;
            }


            if ( B2R_Real_NPatchEachRank[lv][t] == NULL )   B2R_Real_NPatchEachRank[lv][t] = new int [NRank];
            if ( B2R_Buff_NPatchEachRank[lv][t] == NULL )   B2R_Buff_NPatchEachRank[lv][t] = new int [NRank];

            for (int r=0; r<NRank; r++)
            {
               B2R_Real_NPatchEachRank[lv][t][r] = 0;
               B2R_Buff_NPatchEachRank[lv][t][r] = 0;
            }
         } // for (int t=0; t<2; t++)


         if ( F2S_Send_NPatchEachRank[lv] == NULL )   F2S_Send_NPatchEachRank[lv] = new int [NRank];
         if ( F2S_Recv_NPatchEachRank[lv] == NULL )   F2S_Recv_NPatchEachRank[lv] = new int [NRank];

         for (int r=0; r<NRank; r++)
         {
            F2S_Send_NPatchEachRank[lv][r] = 0;
            F2S_Recv_NPatchEachRank[lv][r] = 0;
         }
      } // for (int lv=0; lv<NLEVEL; lv++)
#     endif // #ifdef LOAD_BALANCE

//    set pointers
      Mass = Attribute[PAR_MASS];
      PosX = Attribute[PAR_POSX];
      PosY = Attribute[PAR_POSY];
      PosZ = Attribute[PAR_POSZ];
      VelX = Attribute[PAR_VELX];
      VelY = Attribute[PAR_VELY];
      VelZ = Attribute[PAR_VELZ];
      Time = Attribute[PAR_TIME];
      Type = Attribute[PAR_TYPE];
#     ifdef STORE_PAR_ACC
      AccX = Attribute[PAR_ACCX];
      AccY = Attribute[PAR_ACCY];
      AccZ = Attribute[PAR_ACCZ];
#     endif

   } // METHOD : InitRepo



   //===================================================================================
   // Method      :  AddOneParticle
   // Description :  Add ONE new particle into the particle list
   //
   // Note        :  1. This function will modify several global variables
   //                   --> One must be careful when implementing OpenMP to avoid data race
   //                   --> For example, use the **critical** construct
   //                2. Please provide all necessary information of the new particle
   //                3. If there are inactive particles, their IDs will be reassigned
   //                   to the newly added particles
   //                4. Note that the global variable "AveDensity_Init" will NOT be recalculated
   //                   automatically here
   //
   // Parameter   :  NewAtt : Array storing the attributes of new particles
   //
   // Return      :  Index of the new particle (ParID)
   //===================================================================================
   long AddOneParticle( const real *NewAtt )
   {

//    check
#     ifdef DEBUG_PARTICLE
      if ( NPar_AcPlusInac < 0 ) Aux_Error( ERROR_INFO, "NPar_AcPlusInac (%ld) < 0 !!\n", NPar_AcPlusInac );

      if ( NewAtt == NULL )   Aux_Error( ERROR_INFO, "NewAtt == NULL !!\n" );

      if ( NewAtt[PAR_MASS] < (real)0.0 )
         Aux_Error( ERROR_INFO, "Adding an inactive particle (mass = %21.14e) !!\n", NewAtt[PAR_MASS] );

      if ( NewAtt[PAR_POSX] != NewAtt[PAR_POSX] ||
           NewAtt[PAR_POSY] != NewAtt[PAR_POSY] ||
           NewAtt[PAR_POSZ] != NewAtt[PAR_POSZ]   )
         Aux_Error( ERROR_INFO, "Adding a particle with strange position (%21.14e, %21.14e, %21.14e) !!\n",
                    NewAtt[PAR_POSX], NewAtt[PAR_POSY], NewAtt[PAR_POSZ] );
#     endif


//    1. determine the target particle ID
      long ParID;

//    1-1. reuse an inactive particle ID
      if ( NPar_Inactive > 0 )
      {
         ParID = InactiveParList[ NPar_Inactive-1 ];
         NPar_Inactive --;

#        ifdef DEBUG_PARTICLE
         if ( ParID < 0  ||  ParID >= NPar_AcPlusInac )
            Aux_Error( ERROR_INFO, "Incorrect ParID (%ld), NPar_AcPlusInac = %ld !!\n", ParID, NPar_AcPlusInac );
#        endif
      }

//    1-2. add a new particle ID
      else
      {
//       allocate enough memory for the particle variable array
         if ( NPar_AcPlusInac >= ParListSize )
         {
            ParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*(ParListSize+1) );

            for (int v=0; v<PAR_NATT_TOTAL; v++)   Attribute[v] = (real*)realloc( Attribute[v], ParListSize*sizeof(real) );

            Mass = Attribute[PAR_MASS];
            PosX = Attribute[PAR_POSX];
            PosY = Attribute[PAR_POSY];
            PosZ = Attribute[PAR_POSZ];
            VelX = Attribute[PAR_VELX];
            VelY = Attribute[PAR_VELY];
            VelZ = Attribute[PAR_VELZ];
            Time = Attribute[PAR_TIME];
            Type = Attribute[PAR_TYPE];
#           ifdef STORE_PAR_ACC
            AccX = Attribute[PAR_ACCX];
            AccY = Attribute[PAR_ACCY];
            AccZ = Attribute[PAR_ACCZ];
#           endif
         }

         ParID = NPar_AcPlusInac;
         NPar_AcPlusInac ++;
      } // if ( NPar_Inactive > 0 ) ... else ...


//    2. record the data of new particles
      for (int v=0; v<PAR_NATT_TOTAL; v++)   Attribute[v][ParID] = NewAtt[v];


//    3. update the total number of active particles (assuming all new particles are active)
      NPar_Active ++;


//    4. return the new particle index
      return ParID;

   } // METHOD : AddOneParticle



   //===================================================================================
   // Method      :  RemoveOneParticle
   // Description :  Remove ONE particle from the particle list
   //
   // Note        :  1. This function will modify NPar_Active and NPar_Inactive
   //                   --> Since these are global variables, one must be careful for the OpenMP
   //                       implementation to avoid data race
   //                   --> For example, use the **critical** construct
   //                2. Particles being removed will NOT be actually deleted in the memory.
   //                   Instead, their masses are assigned to Marker, which must be negative
   //                   in order to be distinguishable from active particles.
   //                   --> IDs of all removed particles will be recorded in InactiveParList[]
   //                   --> These IDs will be reassigned to new particles added later on when
   //                       calling AddOneParticle()
   //                3. Note that the global variable "AveDensity_Init" will NOT be recalculated
   //                   automatically here
   //
   // Parameter   :  ParID  : Particle ID to be removed
   //                Marker : Value assigned to the mass of the particle being removed
   //                         (PAR_INACTIVE_OUTSIDE or PAR_INACTIVE_MPI)
   //
   // Return      :  None
   //===================================================================================
   void RemoveOneParticle( const long ParID, const real Marker )
   {

//    check
#     ifdef DEBUG_PARTICLE
      if ( ParID < 0  ||  ParID >= NPar_AcPlusInac )
         Aux_Error( ERROR_INFO, "Wrong ParID (%ld) !!\n", ParID );

      if ( Marker != PAR_INACTIVE_OUTSIDE  &&  Marker != PAR_INACTIVE_MPI )
         Aux_Error( ERROR_INFO, "Unsupported Marker (%14.7e) !!\n", Marker );
#     endif


//    1. allocate enough memory for InactiveParList
      if ( NPar_Inactive >= InactiveParListSize )
      {
         InactiveParListSize = (int)ceil( PARLIST_GROWTH_FACTOR*(InactiveParListSize+1) );

         InactiveParList = (long*)realloc( InactiveParList, InactiveParListSize*sizeof(long) );
      }


//    2. record the particle ID to be removed
      InactiveParList[NPar_Inactive] = ParID;


//    3. remove the target particle
      Mass[ParID] = Marker;

      NPar_Active   --;
      NPar_Inactive ++;

#     ifdef DEBUG_PARTICLE
      if ( NPar_Active + NPar_Inactive != NPar_AcPlusInac )
         Aux_Error( ERROR_INFO, "NPar_Active (%ld) + NPar_Inactive (%ld) != NPar_AcPlusInac (%ld) !!\n",
                    NPar_Active, NPar_Inactive, NPar_AcPlusInac );
#     endif

   } // METHOD : RemoveOneParticle


}; // struct Particle_t



#endif // #ifndef __PARTICLE_H__
