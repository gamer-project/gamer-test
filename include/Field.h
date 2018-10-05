#ifndef __FIELD_H__
#define __FIELD_H__



#include "Macro.h"
#include "Typedef.h"

const FieldIdx_t Idx_Undefined = -1;   // must be outside the range [0, NCOMP_TOTAL-1]


// ***********************************************************
// ** Add new field indices below with the format:          **
// **                                                       **
// ** SET_GLOBAL( FieldIdx_t Idx_NewField, Idx_Undefined ); **
// **                                                       **
// ***********************************************************

#if   ( MODEL == HYDRO )
// main fields
SET_GLOBAL( FieldIdx_t Idx_Dens,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomX,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomY,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomZ,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Engy,          Idx_Undefined );
#ifdef COSMIC_RAY
SET_GLOBAL( FieldIdx_t Idx_CRay,          Idx_Undefined );
#endif
#if   ( DUAL_ENERGY == DE_ENPY )
SET_GLOBAL( FieldIdx_t Idx_Enpy,          Idx_Undefined );
#elif ( DUAL_ENERGY == DE_EINT )
SET_GLOBAL( FieldIdx_t Idx_Eint,          Idx_Undefined );
#endif

// Grackle fields
SET_GLOBAL( FieldIdx_t Idx_e,             Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HI,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HII,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeI,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeII,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeIII,         Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HM,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_H2I,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_H2II,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_DI,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_DII,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HDI,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Metal,         Idx_Undefined );


#elif ( MODEL == ELBDM )

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL



// ***************************************************************
// ** Add new particle attribute indices below with the format: **
// **                                                           **
// ** SET_GLOBAL( FieldIdx_t Idx_NewAttribute, Idx_Undefined ); **
// **                                                           **
// ***************************************************************

#ifdef PARTICLE
SET_GLOBAL( FieldIdx_t Idx_ParMass,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosZ,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelZ,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParType,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParTime,       Idx_Undefined );
#ifdef STORE_PAR_ACC
SET_GLOBAL( FieldIdx_t Idx_ParAccX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParAccY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParAccZ,       Idx_Undefined );
#endif
#ifdef STAR_FORMATION
SET_GLOBAL( FieldIdx_t Idx_ParCreTime,    Idx_Undefined );
#endif
SET_GLOBAL( FieldIdx_t Idx_ParMetalFrac,  Idx_Undefined );
#endif // #ifdef PARTICLE



// field and particle attribute labels
SET_GLOBAL( char FieldLabel[NCOMP_TOTAL][MAX_STRING] );
#ifdef GRAVITY
SET_GLOBAL( char *PotLabel, "Pote" );  // potential label is currently fixed
#endif
#ifdef MHD
SET_GLOBAL( char *MagLabel[NCOMP_MAG] );
#endif
#ifdef PARTICLE
SET_GLOBAL( char ParAttLabel[PAR_NATT_TOTAL][MAX_STRING] );
#endif



#endif  // #ifndef __FIELD_H__
