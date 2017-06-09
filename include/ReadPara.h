#ifndef __READPARA_H__
#define __READPARA_H__


#include <typeinfo>
#include "Global.h"

bool Aux_CheckFileExist( const char *FileName );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );

#define NPARA_MAX    1000     // maximum number of parameters
#define PARA_SET      '#'     // character to indicate that a given parameter has been set properly
#define TYPE_INT        1     // different data types
#define TYPE_LONG       2
#define TYPE_UINT       3
#define TYPE_ULONG      4
#define TYPE_BOOL       5
#define TYPE_FLOAT      6
#define TYPE_DOUBLE     7
#define TYPE_STRING     8




//-------------------------------------------------------------------------------------------------------
// Structure   :  ReadPara_t
// Description :  Data structure ...
//
// Data Member :
//
// Method      :  ReadPara_t : Constructor
//               ~ReadPara_t : Destructor
//                Add        : Add a new parameter
//                Read       : Read parameters from the parameter file
//-------------------------------------------------------------------------------------------------------
struct ReadPara_t
{

// data members
// ===================================================================================
   int    NPara;
   char (*Key)[MAX_STRING];
   void **Ptr;
   int   *Type;
   void **Default;


   //===================================================================================
   // Constructor :  ReadPara_t
   // Description :  Constructor of the structure "ReadPara_t"
   //
   // Note        :  Initialize variables and allocate memory
   //===================================================================================
   ReadPara_t()
   {

      NPara   = 0;
      Key     = new char  [NPARA_MAX][MAX_STRING];
      Ptr     = new void* [NPARA_MAX];
      Type    = new int   [NPARA_MAX];
      Default = new void* [NPARA_MAX];

   } // METHOD : ReadPara_t



   //===================================================================================
   // Constructor :  ~ReadPara_t
   // Description :  Destructor of the structure "ReadPara_t"
   //
   // Note        :  Deallocate memory
   //===================================================================================
   ~ReadPara_t()
   {

      for (int t=0; t<NPara; t++)
      {
         delete Default[t];
      }

      delete [] Key;
      delete [] Ptr;
      delete [] Type;
      delete [] Default;

   } // METHOD : ~ReadPara_t



   //===================================================================================
   // Constructor :  Add
   // Description :  Add a new parameter to be loaded later
   //
   // Note        :  1. This function stores the name, address, and data type of the new parameter
   //                2. Data type (e.g., integer, float, ...) is determined by the input pointer
   //                3. NewPtr, NewDefault, NewMin, and NewMax must have the same data type
   //===================================================================================
   template <typename T>
   void Add( const char NewKey[], T* NewPtr, T NewDefault )
   {

      if ( NPara >= NPARA_MAX )  Aux_Error( ERROR_INFO, "exceed the maximum number of parameters (%d) !!\n", NPARA_MAX );


//    parameter name
      strcpy( Key[NPara], NewKey );

//    parameter address
      Ptr[NPara] = NewPtr;

//    parameter data type
      if      ( typeid(T) == typeid(int   ) )   Type[NPara] = TYPE_INT;
      else if ( typeid(T) == typeid(long  ) )   Type[NPara] = TYPE_LONG;
      else if ( typeid(T) == typeid(uint  ) )   Type[NPara] = TYPE_UINT;
      else if ( typeid(T) == typeid(ulong ) )   Type[NPara] = TYPE_ULONG;
      else if ( typeid(T) == typeid(bool  ) )   Type[NPara] = TYPE_BOOL;
      else if ( typeid(T) == typeid(float ) )   Type[NPara] = TYPE_FLOAT;
      else if ( typeid(T) == typeid(double) )   Type[NPara] = TYPE_DOUBLE;
      else if ( typeid(T) == typeid(char  ) )   Type[NPara] = TYPE_STRING;
      else
         Aux_Error( ERROR_INFO, "unsupported data type for \"%s\" (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n",
                    NewKey );

//    set the default values
      Default[NPara] = new T;
      *( (T*)Default[NPara] ) = NewDefault;

      NPara ++;

   } // METHOD : Add


   //===================================================================================
   // Constructor :  Read
   // Description :  Read all parameters added by Add()
   //
   // Note        :  1. Format:   KEY   VALUE
   //                2. Use # to comment out lines
   //===================================================================================
   void Read( const char *FileName )
   {

      if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );


      char LoadKey[MAX_STRING], LoadValue[MAX_STRING];
      int  MatchIdx, LineNum=0, NLoad;

      char *Line = new char [MAX_STRING];
      FILE *File = fopen( FileName, "r" );


//    loop over all lines in the target file
      while ( fgets( Line, MAX_STRING, File ) != NULL )
      {
         LineNum ++;

//       load the key and value at the target line
         NLoad = sscanf( Line, "%s%s", LoadKey, LoadValue );

//       skip lines with incorrect format (e.g., empyt lines)
         if ( NLoad < 2  &&  MPI_Rank == 0 )
         {
            if ( NLoad == 1 )
               Aux_Error( ERROR_INFO, "cannot find the value to assign to the key \"%s\" at line %d !!\n",
                          LoadKey, LineNum );
            continue;
         }

//       skip lines starting with #
         if ( LoadKey[0] == '#' )   continue;

//       locate the target key
         MatchIdx = -1;
         for (int k=0; k<NPara; k++)
         {
            if (  strcmp( Key[k], LoadKey ) == 0  )
            {
               MatchIdx  = k;
               Key[k][0] = PARA_SET;   // reset the key to any strange character to detect duplicate parameters
               break;
            }
         }

         if ( MatchIdx >= 0 )
         {
            switch ( Type[MatchIdx] )
            {
               case TYPE_INT    :   *( (int*   )Ptr[MatchIdx] ) = (int   )atol( LoadValue );    break;
               case TYPE_LONG   :   *( (long*  )Ptr[MatchIdx] ) = (long  )atol( LoadValue );    break;
               case TYPE_UINT   :   *( (uint*  )Ptr[MatchIdx] ) = (uint  )atol( LoadValue );    break;
               case TYPE_ULONG  :   *( (ulong* )Ptr[MatchIdx] ) = (ulong )atol( LoadValue );    break;
               case TYPE_BOOL   :   *( (bool*  )Ptr[MatchIdx] ) = (bool  )atol( LoadValue );    break;
               case TYPE_FLOAT  :   *( (float* )Ptr[MatchIdx] ) = (float )atof( LoadValue );    break;
               case TYPE_DOUBLE :   *( (double*)Ptr[MatchIdx] ) = (double)atof( LoadValue );    break;
               case TYPE_STRING :   strcpy( (char*)Ptr[MatchIdx], LoadValue );                  break;

               default: Aux_Error( ERROR_INFO, "unsupported data type (char*, float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n" );
            }
         }

         else
         {
            Aux_Error( ERROR_INFO, "unrecognizable parameter \"%s\" at line %d (either non-existing or duplicate) !!\n",
                       LoadKey, LineNum );
         }
      } // while ( fgets( Line, MAX_STRING, File ) != NULL )


      fclose( File );
      delete [] Line;


//    set the parameters missing in the runtime parameter file to their default values
      SetDefault();


//    validate the parameters


   } // METHOD : Read


   //===================================================================================
   // Constructor :  SetDefault
   // Description :  Set parameters missing in the runtime parameter file to their default values
   //
   // Note        :  1. Parameters already set by the runtime parameter file have Key[0] == PARA_SET
   //                2. We do NOT set default values for strings. An error will be raised instead.
   //===================================================================================
   void SetDefault()
   {

      for (int t=0; t<NPara; t++)
      {
         if ( Key[t][0] != PARA_SET )
         {
            long   def_int = NULL_INT;
            double def_flt = NULL_REAL;

            switch ( Type[t] )
            {
               case TYPE_INT    :   def_int = *( (int*   )Ptr[t] ) = *( (int*   )Default[t] );  break;
               case TYPE_LONG   :   def_int = *( (long*  )Ptr[t] ) = *( (long*  )Default[t] );  break;
               case TYPE_UINT   :   def_int = *( (uint*  )Ptr[t] ) = *( (uint*  )Default[t] );  break;
               case TYPE_ULONG  :   def_int = *( (ulong* )Ptr[t] ) = *( (ulong* )Default[t] );  break;
               case TYPE_BOOL   :   def_int = *( (bool*  )Ptr[t] ) = *( (bool*  )Default[t] );  break;
               case TYPE_FLOAT  :   def_flt = *( (float* )Ptr[t] ) = *( (float* )Default[t] );  break;
               case TYPE_DOUBLE :   def_flt = *( (double*)Ptr[t] ) = *( (double*)Default[t] );  break;
               case TYPE_STRING :   Aux_Error( ERROR_INFO, "string variable \"%s\" is not set !!\n", Key[t] ); break;

               default: Aux_Error( ERROR_INFO, "no default value for this data type (float*, double*, int*, long*, unit*, ulong*, bool* only) !!\n" );
            }

            if ( MPI_Rank == 0 )
            {

               if ( def_int != NULL_INT )
                  Aux_Message( stdout, "NOTE : parameter [%20s] is set to the default value [%21ld]\n",
                               Key[t], def_int );
               else
                  Aux_Message( stdout, "NOTE : parameter [%20s] is set to the default value [%21.14e]\n",
                               Key[t], def_flt );
            }
         }
      } // for (int t=0; t<NPara; t++)

   } // METHOD : SetDefault

}; // struct ReadPara_t



// remove symbolic constants since they are only used in this structure
#undef NPARA_MAX
#undef PARA_SET
#undef TYPE_INT
#undef TYPE_LONG
#undef TYPE_UINT
#undef TYPE_ULONG
#undef TYPE_FLOAT
#undef TYPE_DOUBLE
#undef TYPE_BOOL
#undef TYPE_STRING



#endif // #ifndef __READPARA_H__
