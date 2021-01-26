#ifndef __PHYSICAL_CONSTANT_H__
#define __PHYSICAL_CONSTANT_H__



// *********************************************************************************
// ** This header defines common physical constants.                              **
// ** All values are in CGS units unless otherwise specified.                     **
// **                                                                             **
// ** Reference: K.A. Olive et al. (Particle Data Group),                         **
// ** Chin. Phys. C, 38, 090001 (2014) and 2015 update.                           **
// ** (1) http://pdg.lbl.gov/2015/reviews/rpp2015-rev-phys-constants.pdf          **
// ** (2) http://pdg.lbl.gov/2015/reviews/rpp2015-rev-astrophysical-constants.pdf **
// *********************************************************************************


// length units
const double Const_cm            = 1.0;
const double Const_m             = 1.0e2*Const_cm;
const double Const_km            = 1.0e5*Const_cm;
const double Const_pc            = 3.08567758149e18;        // parsec
const double Const_kpc           = 1.0e3*Const_pc;
const double Const_Mpc           = 1.0e6*Const_pc;
const double Const_Gpc           = 1.0e9*Const_pc;
const double Const_au            = 1.495978707e13;          // astronomical unit
const double Const_ly            = 9.46053e17;              // light year

// mass units
const double Const_g             = 1.0;
const double Const_kg            = 1.0e3*Const_g;
const double Const_Msun          = 1.9885e33;               // solar mass
const double Const_amu           = 1.660539040e-24;         // atomic mass unit
const double Const_mH            = 1.007825*Const_amu;      // neutral atomic hydrogen mass
const double Const_mp            = 1.672621898e-24;         // proton mass
const double Const_me            = 9.10938356e-28;          // electron mass

// time units
const double Const_s             = 1.0;                     // second
const double Const_yr            = 3.15569252e7;            // year
const double Const_kyr           = 1.0e3*Const_yr;
const double Const_Myr           = 1.0e6*Const_yr;
const double Const_Gyr           = 1.0e9*Const_yr;

// temperature units
const double Const_K             = 1.0;                     // kelvin

// magnetic field units
const double Const_Gauss         = 1.0;                     // gauss
const double Const_mGauss        = 1.0e-3;                  // milligauss
const double Const_muGauss       = 1.0e-6;                  // microgauss
const double Const_nGauss        = 1.0e-9;                  // nanogauss
const double Const_Tesla         = 1.0e4;                   // tesla
const double Const_mTesla        = 1.0e1;                   // millitesla
const double Const_muTesla       = 1.0e-2;                  // microtesla
const double Const_nTesla        = 1.0e-5;                  // nanotesla

// other physical constants
const double Const_c             = 2.99792458e10;           // speed of light
const double Const_Planck        = 1.054571800e-27;         // reduced Planck constant in erg*s
const double Const_Planck_eV     = 6.582119514e-16;         // reduced Planck constant in eV*s
const double Const_NewtonG       = 6.6738e-8;               // gravitational constant
const double Const_eV            = 1.6021766208e-12;        // electron volt
const double Const_keV           = 1.0e3*Const_eV;
const double Const_MeV           = 1.0e6*Const_eV;
const double Const_GeV           = 1.0e9*Const_eV;
const double Const_kB            = 1.38064852e-16;          // Boltzmann constant in erg/K
const double Const_kB_eV         = 8.6173303e-5;            // Boltzmann constant in eV/K
const double Const_NA            = 6.022140857e23;          // Avogadro constant (in 1/mole)



#endif  // #ifndef __PHYSICAL_CONSTANT_H__
