#include "GAMER.h"

#ifdef PARTICLE


// declare as static so that other functions cannot invoke it directly and must use the function pointer
void Par_Init_Attribute_User_Template();

// this function pointer must be set by a test problem initializer
void (*Par_Init_Attribute_User_Ptr)() = NULL;


static int NDefinedAtt;    // total number of defined attributes --> for debug only




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute
// Description :  Set all particle attributes (including both built-in and user-defined attributes)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Total number of attributes is determined by PAR_NATT_TOTAL = PAR_NATT_BUILTIN + PAR_NATT_USER
//                3. To initialize user-defined particle attributes, the function pointer "Par_Init_Attribute_User_Ptr"
//                   must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 0. initialize counter as zero
   NDefinedAtt = 0;


// 1. add built-in intrinsic attributes
//    --> must not change the following order of declaration since they must be consistent
//        with the symbolic constants defined in Macro.h (e.g., PAR_MASS)
   Idx_ParMass    = AddParticleAttribute( "ParMass" );
   Idx_ParPosX    = AddParticleAttribute( "ParPosX" );
   Idx_ParPosY    = AddParticleAttribute( "ParPosY" );
   Idx_ParPosZ    = AddParticleAttribute( "ParPosZ" );
   Idx_ParVelX    = AddParticleAttribute( "ParVelX" );
   Idx_ParVelY    = AddParticleAttribute( "ParVelY" );
   Idx_ParVelZ    = AddParticleAttribute( "ParVelZ" );
   Idx_ParType    = AddParticleAttribute( "ParType" );


// 2. add other built-in attributes
#  ifdef STAR_FORMATION
   Idx_ParCreTime = AddParticleAttribute( "ParCreTime" );
#  endif


// 3. add user-defined attributes
   if ( Par_Init_Attribute_User_Ptr != NULL )   Par_Init_Attribute_User_Ptr();


// 4. must put built-in attributes not to be stored on disk at the END of the attribute list
//    --> make it easier to discard them when storing data on disk (see Output_DumpData_Total(_HDF5).cpp)
//    --> must also be consistent with the symbolic constant (e.g., PAR_TIME and PAR_ACC*) defined in Macro.h
//    --> total number of attributes not to be stored on disk is set by PAR_NATT_UNSTORED
//        --> currently including time and acceleration*3
#  ifdef STORE_PAR_ACC
   Idx_ParAccX    = AddParticleAttribute( "ParAccX" );
   Idx_ParAccY    = AddParticleAttribute( "ParAccY" );
   Idx_ParAccZ    = AddParticleAttribute( "ParAccZ" );
#  endif
   Idx_ParTime    = AddParticleAttribute( "ParTime" );



// 5. validate if all attributes have been set properly
   if ( NDefinedAtt != PAR_NATT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined attributes (%d) != expectation (%d) !!\n"
                 "        --> Modify PAR_NATT_USER in the Makefile or invoke AddParticleAttribute() properly\n",
                 NDefinedAtt, PAR_NATT_TOTAL );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_Attribute



//-------------------------------------------------------------------------------------------------------
// Function    :  AddParticleAttribute
// Description :  Add a new particle attribute to the attribute list
//
// Note        :  1. This function will
//                   (1) set the attribute label, which will be used as the output name of the attribute
//                   (2) return the index of the new attribute, which can be used to access the attribute
//                       data (e.g., amr->Par->Attribute[])
//                2. One must invoke AddParticleAttribute() exactly PAR_NATT_TOTAL times to set the labels
//                   of all attributes
//                3. Invoked by Par_Init_Attribute() and various test problem initializers
//
// Parameter   :  InputLabel : Label (i.e., name) of the new attribute
//
// Return      :  (1) ParAttLabel[]
//                (2) Index of the newly added attribute
//-------------------------------------------------------------------------------------------------------
FieldIdx_t AddParticleAttribute( char *InputLabel )
{

   const FieldIdx_t AttIdx = NDefinedAtt ++;


// check
   if ( InputLabel == NULL )
      Aux_Error( ERROR_INFO, "InputLabel == NULL !!\n" );

   if ( NDefinedAtt > PAR_NATT_TOTAL )
      Aux_Error( ERROR_INFO, "total number of defined particle attributes (%d) exceeds expectation (%d) after adding \"%s\" !!\n"
                 "        --> Modify PAR_NATT_USER in the Makefile properly\n",
                 NDefinedAtt, PAR_NATT_TOTAL, InputLabel );

   for (int v=0; v<NDefinedAtt-1; v++)
      if (  strcmp( ParAttLabel[v], InputLabel ) == 0  )
         Aux_Error( ERROR_INFO, "duplicate particle attribute label \"%s\" !!\n", InputLabel );


// set attribute label
   strcpy( ParAttLabel[AttIdx], InputLabel );


// return attribute index
   return AttIdx;

} // FUNCTION : AddParticleAttribute



//-------------------------------------------------------------------------------------------------------
// Function    :  GetParticleAttributeIndex
// Description :  Return the index of the target particle attribute
//
// Note        :  1. Usage: AttIdx = GetParticleAttributeIndex( ParAttLabel );
//                2. Return Idx_Undefined if the target attribute cannot be found
//
// Parameter   :  InputLabel : Target attribute label
//                Check      : Whether or not to terminate the program if the target attribute cannot be found
//                             --> Accepted options: CHECK_ON / CHECK_OFF
//
// Return      :  Sucess: index of the target attribute
//                Failed: Idx_Undefined
//-------------------------------------------------------------------------------------------------------
FieldIdx_t GetParticleAttributeIndex( char *InputLabel, const Check_t Check )
{

   FieldIdx_t Idx_Out = Idx_Undefined;

   for (int v=0; v<NDefinedAtt; v++)
   {
      if (  strcmp( ParAttLabel[v], InputLabel ) == 0  )
      {
         Idx_Out = v;
         break;
      }
   }

   if ( Check == CHECK_ON  &&  Idx_Out == Idx_Undefined )
      Aux_Error( ERROR_INFO, "cannot find the target particle attribute \"%s\" !!\n", InputLabel );

   return Idx_Out;

} // FUNCTION : GetParticleAttributeIndex



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_Attribute_User_Template
// Description :  Template of adding user-defined particle attributes
//
// Note        :  1. Invoked by Par_Init_Attribute() using the function pointer "Par_Init_Attribute_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Init_Attribute_User_Template()
{

// example
// Idx_NewAtt = AddParticleAttribute( "NewAttLabel" );

} // FUNCTION : Par_Init_Attribute_User_Template



#endif // #ifdef PARTICLE
