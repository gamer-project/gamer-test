#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Scale2PhySize
// Description :  Scale --> corresponding physical size
//
// Parameter   :  Scale :  Input scale
//
// Return      :  Corresponding physical size
//-------------------------------------------------------------------------------------------------------
double Mis_Scale2PhySize( const int Scale )
{

   return Scale*amr->dh[TOP_LEVEL];

} // FUNCTION : Mis_Scale2PhySize



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Cell2PhySize
// Description :  Number of cells at the target level --> corresponding physical size
//
// Parameter   :  NCell :  Input number of cells
//                lv    :  Target level
//
// Return      :  Corresponding physical size
//-------------------------------------------------------------------------------------------------------
double Mis_Cell2PhySize( const int NCell, const int lv )
{

   return NCell*amr->dh[lv];

} // FUNCTION : Mis_Cell2PhySize



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Scale2Cell
// Description :  Scale --> corresponding number of cells at the target level
//
// Parameter   :  Scale :  Input scale
//                lv    :  Target level
//
// Return      :  Corresponding number of cells at lv
//-------------------------------------------------------------------------------------------------------
int Mis_Scale2Cell( const int Scale, const int lv )
{

#  ifdef GAMER_DEBUG
   if ( Scale % amr->scale[lv] != 0 )
      Aux_Message( stderr, "WARNING : input scale (%d) is not a multiple of the scale unit (lv %d --> %d) !!\n",
                   Scale, lv, amr->scale[lv] );
#  endif

   return ( Scale >> (TOP_LEVEL-lv) );

} // FUNCTION : Mis_Scale2Cell



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Cell2Scale
// Description :  Number of cells at the target level --> scale
//
// Parameter   :  NCell :  Input number of cells
//                lv    :  Target level
//
// Return      :  Corresponding scale
//-------------------------------------------------------------------------------------------------------
int Mis_Cell2Scale( const int NCell, const int lv )
{

   return NCell*amr->scale[lv];

} // FUNCTION : Mis_Cell2Scale



//-------------------------------------------------------------------------------------------------------
// Function    : Mis_Cartesian2Spherical
// Description : Convert Cartesian coordinates to spherical coordinates
//
// Note        : Reference: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
//
// Parameter   : Cartesian    : The array stored Cartesian coordinates
//                              --> Cartesian[0/1/2] = x/y/z
//               Spherical    : The array stored spherical coordinates
//                              --> Spherical[0] = r, the length of radial vector
//                              --> Spherical[1] = θ, the angle between the z-axis and the radial vector
//                              --> Spherical[2] = φ, the angle between the x-axis and the projection of
//                                                    the radial vector onto the xy-plane
// Return      : Spherical[]
//-------------------------------------------------------------------------------------------------------
void Mis_Cartesian2Spherical( const double Cartesian[], double Spherical[] )
{

  Spherical[0] = SQRT( SQR(Cartesian[0]) + SQR(Cartesian[1]) + SQR(Cartesian[2]) );

  if ( SQR(Cartesian[0]) + SQR(Cartesian[1]) == 0.0 && Cartesian[2] == 0.0 )
    Aux_Error( ERROR_INFO, "Both arguments in atan2 can not be zero !! (%s)\n", __FUNCTION__ );

  Spherical[1] = ATAN2( SQRT( SQR(Cartesian[0]) + SQR(Cartesian[1]) ), Cartesian[2] );

  if ( Cartesian[0] == 0.0 && Cartesian[1] == 0.0 )
    Aux_Error( ERROR_INFO, "Both arguments in atan2 can not be zero !! (%s)\n", __FUNCTION__ );

  Spherical[2] = ATAN2( Cartesian[1], Cartesian[0] );

}



//-------------------------------------------------------------------------------------------------------
// Function    : Mis_Cartesian2Cylindrical
// Description : Convert Cartesian coordinates to cylindrical coordinates
//
// Note        : Reference: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
//
// Parameter   : Cartesian    : The array stored Cartesian coordinates
//                              --> Cartesian[0/1/2] = x/y/z
//               Cylindrical  : The array stored cylindrical coordinates
//                            --> Cylindrical[0] = ρ, the length of radial vector
//                            --> Cylindrical[1] = φ, the angle between the x-axis and the radial vector
//                            --> Cylindrical[2] = z, z coordinates
// Return      : Cylindrical[]
//-------------------------------------------------------------------------------------------------------
void Mis_Cartesian2Cylindrical( const double Cartesian[], double Cylindrical[] )
{

  Cylindrical[0] = SQRT( SQR(Cartesian[0]) + SQR(Cartesian[1]) );

  if ( Cartesian[1] == 0.0 && Cartesian[0] == 0.0 )
    Aux_Error( ERROR_INFO, "Both arguments in atan2 can not be zero !! (%s)\n", __FUNCTION__ );

  Cylindrical[1] = ATAN2( Cartesian[1], Cartesian[0] );

  Cylindrical[2] = Cartesian[2];

}



//-------------------------------------------------------------------------------------------------------
// Function    : Mis_Spherical2Cartesian
// Description : Convert Spherical coordinates to Cartesian coordinates
//
// Note        : Reference: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
//
//
// Parameter   : Spherical    : The array stored Spherical coordinates
//                              --> Spherical[0] = r, the length of radial vector
//                              --> Spherical[1] = θ, the angle between the z-axis and the radial vector
//                              --> Spherical[2] = φ, the angle between the x-axis and the projection of
//                                                    the radial vector onto the xy-plane
//               Cartesian    : The array stored Cartesian coordinates
//                              --> Cartesian[0/1/2] = x/y/z
//
// Return      : Cartesian[]
//-------------------------------------------------------------------------------------------------------
void Mis_Spherical2Cartesian( const double Spherical[], double Cartesian[] )
{

  Cartesian[0] = Spherical[0]*SIN(Spherical[1])*COS(Spherical[2]);

  Cartesian[1] = Spherical[0]*SIN(Spherical[1])*SIN(Spherical[2]);

  Cartesian[2] = Spherical[0]*COS(Spherical[1]);

}



//-------------------------------------------------------------------------------------------------------
// Function    : Mis_RotateRigidBody
// Description : 1. Rotate rigid body based on Euler angles.
//               2. Rotation convention follows the section 4.4 in Classical Mechenial, Goldstein, 3rd edition
//                  --> Using Eq.(4.47)
//               3. Rotating rigid body is equivalent to active rotation
//               4. However, below process is for passive coordinate transformation
//               5. First, rotating the initial axes, xyz, by an angle φ counterclockwise about the z-axis.
//                  The resultant coordinate system is labeled the ξηζ axes.
//                  Second, rotate ξηζ axes about the ξ-axis counterclockwise by an angle θ to produce another
//                  intermediate set, ξ'η'ζ' axes.
//                  Finally, the ξ'η'ζ' axes are rotated counterclockwise by an angle ψ about the ζ' axis to produce
//                  the desired x'y'z' system.
//               6. The angle set (φ, θ, ψ) thus completely specify the 3D rotation from the xyz to x'y'z' system.
//
// Note        : 1. This function rotate and change the array `Cartesian` in-place
//               2. The Euler angles between two given systems can be found by ....
//
// Parameter   : Cartesian         : The array stored Cartesian coordinates
//                                  --> Cartesian[0/1/2] = x/y/z
//               EulerAngle[0/1/2] : φ/θ/ψ
//
//
// Return      : Cartesian[]
//-------------------------------------------------------------------------------------------------------
void Mis_RotateRigidBody( double Cartesian[], const double EulerAngle[] )
{

  double CartesianRot[3], Phi, Theta, Psi;

  Phi   = EulerAngle[0];
  Theta = EulerAngle[1];
  Psi   = EulerAngle[2];

  CartesianRot[0] =    ( COS(Phi)*COS(Psi) - COS(Theta)*SIN(Phi)*SIN(Psi) )*Cartesian[0]
                     - ( COS(Phi)*SIN(Psi) + COS(Theta)*SIN(Phi)*COS(Psi) )*Cartesian[1]
                     + (                     SIN(Theta)*SIN(Phi)          )*Cartesian[2];

  CartesianRot[1] =  - ( SIN(Phi)*COS(Psi) + COS(Theta)*COS(Phi)*SIN(Psi) )*Cartesian[0]
                     - ( SIN(Phi)*SIN(Psi) - COS(Theta)*COS(Phi)*COS(Psi) )*Cartesian[1]
                     - (                     SIN(Theta)*COS(Phi)          )*Cartesian[2];

  CartesianRot[2] = SIN(Theta)*SIN(Psi)*Cartesian[0] + SIN(Theta)*COS(Psi)*Cartesian[1] + COS(Theta)*Cartesian[2];

  for (int idx=0; idx<3; idx++) Cartesian[idx] = CartesianRot[idx];

}



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetEulerAngles
// Description :  Get Euler angles by 
// Note        :  See https://en.wikipedia.org/wiki/Euler_angles#Proper_Euler_angles_2
//
// Parameter   :  EulerAngle[0/1/2] : φ/θ/ψ, Euler angles, described in the function `RotateRigidBody`
//             :  Unit[0]           : The z-coordinate of the unit vector along y'-axis in the original xyz system
//             :  Unit[1]           : The y-coordinate of the unit vector along z'-axis in the original xyz system
//             :  Unit[2]           : The z-coordinate of the unit vector along z'-axis in the original xyz system
//
// Return      :  EulerAngle[0/1/2] : φ/θ/ψ
//-------------------------------------------------------------------------------------------------------

void Mis_GetEulerAngles( double EulerAngle[], const double Unit[] )
{

  double Argument0 = -Unit[1] / SQRT( (double)1.0 - Unit[2]*Unit[2] );
  double Argument1 = +Unit[2];
  double Argument2 = +Unit[0] / SQRT( (double)1.0 - Unit[2]*Unit[2] );

  bool Fail = false;

  Fail |= Argument0 < (double)-1.0 || Argument0 > (double)1.0;
  Fail |= Argument1 < (double)-1.0 || Argument1 > (double)1.0;
  Fail |= Argument2 < (double)-1.0 || Argument2 > (double)1.0;

  if ( Fail )
    Aux_Error( ERROR_INFO, "Argument should be lie between -1 and 1 !! (%s)\n", Argument2, __FUNCTION__ );

  EulerAngle[0] =  acos( Argument0 );
  EulerAngle[1] =  acos( Argument1 );
  EulerAngle[2] =  acos( Argument2 );

}
