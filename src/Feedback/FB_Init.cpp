#include "GAMER.h"

#ifdef FEEDBACK



// prototypes of built-in feedbacks
void FB_Init_SNe();


// user-specified feedback to be set by a test problem initializer
void (*FB_Init_User_Ptr)() = NULL;

extern void (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                            real (*Fluid)[PS2][PS2][PS2], const double EdgeL[], const double dh, bool CoarseFine[],
                            const int TID, RandomNumber_t *RNG );


// random number generators
RandomNumber_t *FB_RNG = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init
// Description :  Initialize feedback
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set "FB_Init_User_Ptr" in a test problem initializer for a user-specified feedback
//                3. Allocate NT thread-safe random number generators, where NT is the number of OpenMP threads
//                   --> Freed by FB_End()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_Init()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check if at least one feedback is activated
   FB_Any = ( FB_SNE || FB_USER );

   if ( ! FB_Any  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : all feedback options (e.g., FB_SNE) are inactivated for FEEDBACK !!\n" );


// call the initialization routines of different feedbacks
   if ( FB_SNE)   FB_Init_SNe();

   if ( FB_USER )
   {
      if ( FB_Init_User_Ptr != NULL )  FB_Init_User_Ptr();
      else                             Aux_Error( ERROR_INFO, "FB_Init_User_Ptr == NULL for FB_USER !!\n" );

//    check if FB_User_Ptr is set properly by FB_Init_User_Ptr()
      if ( FB_User_Ptr == NULL )    Aux_Error( ERROR_INFO, "FB_User_Ptr == NULL for FB_USER !!\n" );
   }


// allocate the random number generators for all OpenMP threads
   int NT;
#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {  NT = omp_get_num_threads();  }
#  else
   {  NT = 1;                      }
#  endif

   FB_RNG = new RandomNumber_t( NT );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : FB_Init



#endif // #ifdef FEEDBACK
