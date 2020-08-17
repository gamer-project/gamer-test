#include "GAMER.h"

#if ( MODEL == HYDRO )



// prototypes of built-in EoS
#if   ( EOS == EOS_GAMMA )
void EoS_Init_Gamma();
#elif ( EOS == EOS_NUCLEAR )
# error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#endif // # EOS




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Init
// Description :  Initialize the equation of state (EoS)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. For a non-built-in EoS, "EoS_Init_Ptr" must be set by a test problem initializer
//                   in advance
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void EoS_Init()
{

// set the initialization function pointer for the built-in EoS
#  if   ( EOS == EOS_GAMMA )
   EoS_Init_Ptr = EoS_Init_Gamma;
#  elif ( EOS == EOS_NUCLEAR )
#  error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#  endif // # EOS


// initialize EoS
   if ( EoS_Init_Ptr != NULL )
      EoS_Init_Ptr();
   else
      Aux_Error( ERROR_INFO, "EoS_Init_Ptr == NULL for EoS %d !!\n", EOS );

} // FUNCTION : EoS_Init



#endif // #if ( MODEL == HYDRO )