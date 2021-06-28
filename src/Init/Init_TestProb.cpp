#include "GAMER.h"


// ******************************************************************
// add the new test problem function prototypes here
// ******************************************************************
void Init_TestProb_Hydro_BlastWave();
void Init_TestProb_Hydro_AcousticWave();
void Init_TestProb_Hydro_Bondi();
void Init_TestProb_Hydro_ClusterMerger();
void Init_TestProb_Hydro_AGORA_IsolatedGalaxy();
void Init_TestProb_Hydro_Caustic();
void Init_TestProb_Hydro_SphericalCollapse();
void Init_TestProb_Hydro_KelvinHelmholtzInstability();
void Init_TestProb_Hydro_Riemann();
void Init_TestProb_Hydro_CollidingJets();
void Init_TestProb_Hydro_Plummer();
void Init_TestProb_Hydro_Gravity();
void Init_TestProb_Hydro_MHD_ABC();
void Init_TestProb_Hydro_MHD_OrszagTangVortex();
void Init_TestProb_Hydro_MHD_LinearWave();
void Init_TestProb_Hydro_JeansInstability();
void Init_TestProb_Hydro_NFW();
void Init_TestProb_Hydro_BarredPot();

void Init_TestProb_ELBDM_ExtPot();




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize the target test problem
//
// Note        :  1. Use TESTPROB_ID to choose the target test problem
//                2. All test problem IDs are defined in "include/Typedef.h"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

// ******************************************************************
// add the new test problem IDs here
// ******************************************************************
   switch ( TESTPROB_ID )
   {
      case TESTPROB_NONE :                                                                                  break;

      case TESTPROB_HYDRO_BLAST_WAVE :                   Init_TestProb_Hydro_BlastWave();                   break;
      case TESTPROB_HYDRO_ACOUSTIC_WAVE :                Init_TestProb_Hydro_AcousticWave();                break;
      case TESTPROB_HYDRO_BONDI :                        Init_TestProb_Hydro_Bondi();                       break;
      case TESTPROB_HYDRO_CLUSTER_MERGER :               Init_TestProb_Hydro_ClusterMerger();               break;
      case TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY :        Init_TestProb_Hydro_AGORA_IsolatedGalaxy();        break;
      case TESTPROB_HYDRO_CAUSTIC :                      Init_TestProb_Hydro_Caustic();                     break;
      case TESTPROB_HYDRO_SPHERICAL_COLLAPSE :           Init_TestProb_Hydro_SphericalCollapse();           break;
      case TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY : Init_TestProb_Hydro_KelvinHelmholtzInstability();  break;
      case TESTPROB_HYDRO_RIEMANN :                      Init_TestProb_Hydro_Riemann();                     break;
      case TESTPROB_HYDRO_COLLIDING_JETS :               Init_TestProb_Hydro_CollidingJets();               break;
      case TESTPROB_HYDRO_PLUMMER :                      Init_TestProb_Hydro_Plummer();                     break;
      case TESTPROB_HYDRO_BARRED_POT :                   Init_TestProb_Hydro_BarredPot();                   break;
      case TESTPROB_HYDRO_GRAVITY :                      Init_TestProb_Hydro_Gravity();                     break;
      case TESTPROB_HYDRO_MHD_ABC :                      Init_TestProb_Hydro_MHD_ABC();                     break;
      case TESTPROB_HYDRO_MHD_ORSZAG_TANG_VORTEX :       Init_TestProb_Hydro_MHD_OrszagTangVortex();        break;
      case TESTPROB_HYDRO_MHD_LINEAR_WAVE :              Init_TestProb_Hydro_MHD_LinearWave();              break;
      case TESTPROB_HYDRO_JEANS_INSTABILITY :            Init_TestProb_Hydro_JeansInstability();            break;
      case TESTPROB_HYDRO_NFW :                          Init_TestProb_Hydro_NFW();                         break;

      case TESTPROB_ELBDM_EXTPOT :                       Init_TestProb_ELBDM_ExtPot();                      break;

      default: Aux_Error( ERROR_INFO, "unsupported TESTPROB_ID (%d) !!\n", TESTPROB_ID );
   } // switch( TESTPROB_ID )

} // FUNCTION : Init_TestProb
