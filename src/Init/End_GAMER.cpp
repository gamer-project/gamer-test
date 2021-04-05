#include "GAMER.h"

extern void (*End_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  End_GAMER
// Description :  Put everything you want to do before terminating the program right here
//
// Note        :  1. Function pointer "End_User_Ptr" may be set by a test problem initializer
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_GAMER()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  ifdef TIMING
   Aux_DeleteTimer();
#  endif

   End_MemFree();

   if ( End_User_Ptr != NULL )   End_User_Ptr();

#  ifdef GRAVITY
   End_FFTW();
#  endif

#  ifdef SUPPORT_LIBYT
   YT_End();
#  endif

#  ifdef SUPPORT_GRACKLE
   Grackle_End();
#  endif

#  if ( MODEL == HYDRO )
   EoS_End();
#  endif

   Src_End();

#  ifdef FEEDBACK
   FB_End();
#  endif

#  ifdef GRAVITY
   End_ExtAccPot();
#  endif


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
      Aux_Message( stdout, "\n\n~ GAME OVER ~\n\n\n" );
   }

   MPI_Finalize();

   exit( 0 );

} // FUNCTION : End_GAMER
