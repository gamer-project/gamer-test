#include <cstdio>
#include <cmath>

typedef double real; // this must be consistent with the compilation option FLOAT8

int main()
{
// parameters
   const char  *Filename = "ExtPotTable";       // EXT_POT_TABLE_NAME
   const int    NX       = 131;                 // EXT_POT_TABLE_NPOINT_X
   const int    NY       = 131;                 // EXT_POT_TABLE_NPOINT_Y
   const int    NZ       = 131;                 // EXT_POT_TABLE_NPOINT_Z
   const double dh       = 0.0078125;           // EXT_POT_TABLE_NPOINT_DH
   const double EdgeL[3] = { -dh, -dh, -dh };   // EXT_POT_TABLE_EDGEL_X/Y/Z

// set the external potential array
   const double G  = 1.0;                       // gravitational constant
   const double M  = 1.0;                       // particle mass
   const double cx = EdgeL[0] + 0.5*(NX-1)*dh;  // centre
   const double cy = EdgeL[1] + 0.5*(NY-1)*dh;
   const double cz = EdgeL[2] + 0.5*(NZ-1)*dh;

   double dx, dy, dz, r;

   real (*ExtPot)[NY][NX] = new real [NZ][NY][NX];

   for (int k=0; k<NZ; k++)  {  dz = EdgeL[2] + k*dh - cz;
   for (int j=0; j<NY; j++)  {  dy = EdgeL[1] + j*dh - cy;
   for (int i=0; i<NX; i++)  {  dx = EdgeL[0] + i*dh - cx;

      r               = sqrt( dx*dx + dy*dy + dz*dz );
      ExtPot[k][j][i] = -G*M/r;

   }}}

// output
   FILE *File = fopen( Filename, "wb" );
   fwrite( ExtPot, sizeof(real), (long)NX*NY*NZ, File );
   fclose( File );

   delete [] ExtPot;
}
