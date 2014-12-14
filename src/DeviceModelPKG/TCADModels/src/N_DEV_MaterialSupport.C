//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MaterialSupport.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/19/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.33 $
//
// Revision Date  : $Date: 2014/04/24 21:57:47 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_MaterialSupport.h>

namespace Xyce {
namespace Device {

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassN
// Purpose       : returns effective mass for electrons.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassN (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double mass=0.0;

  if (mater == "si")
  {
    double ml = 0.98; // longitudinal mass
    double mt = 0.19; // transverse mass
    mass = pow((ml*mt*mt),1.0/3.0);
  }
  else if (mater == "ge" )
  {
    double ml = 1.64; // longitudinal mass
    double mt = 0.082; // transverse mass
    mass = pow((1.64*0.082*0.082),1.0/3.0);
  }
  else if (mater == "gaas")
  {
    //mass = 0.067;
    mass = 6.69935094e-2;
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    mass = 0.074;
    //dnco=  0.020;
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    mass = 0.041;
    // dnco=  0.0083;
  }
  else if (mater == "ingap")
  {
    mass = 0.0179;
    // dnco =  0.02391;
  }
  else if (mater == "inp")
  {
    mass = 0.079;
    // dnco=  0.01734;
  }
  else
  {
    Report::UserFatal0() << material << " material not recognized.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassP
// Purpose       : returns effective mass for holes.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassP (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass=0.0;

  if (mater == "si")
  {
    double mlh = 0.16; // light hole mass
    double mhh = 0.49; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh,1.5)),2.0/3.0);
  }
  else if (mater == "ge" )
  {
    double mlh = 0.04; // light hole mass
    double mhh = 0.28; // heavy hole mass
    mass = pow((pow(mlh, 1.5) + pow(mhh, 1.5)),2.0/3.0);
  }
  else if (mater == "gaas")
  {
#if 1
    // Note: Wampler's value of dnva = 0.4318
    // Their SAND report uses a DOS effective mass as: mh*=0.571
    // Note that the electron effective mass for GaAs seems to be uncontroversial m=0.067
    // In order to match XPD, and pass historical regression tests, etc, using this number.
    //mass=0.571; 
    mass=5.71287987e-1;
#endif
#if 0
    // These are the original values in this file:
    // Streetman, Appendix III:
    //  This gives dnva = 0.373684
    double mlh = 0.074; // light hole mass
    double mhh = 0.5;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
#endif
#if 0
    // Websites:  http://ecee.colorado.edu/~bart/book/book/chapter2/pdf/ch2_3_7.pdf
    // http://ecee.colorado.edu/~bart/book/effmass.htm#dosmass
    // http://www.semiconductors.co.uk/propiiiv5653.htm
    //
    // These have been referenced as coming from "Singh, 1993":  
    // Singh J, 1993, "Physics of Semiconductors and Their Heterostructures" (McGraw-Hill) 
    //
    // This gives dnva =  0.32535
    double mlh = 0.082; // light hole mass
    double mhh = 0.45;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
#endif
#if 0
    // William Liu HBT book:
    // this gives dnva = 0.51385
    double mlh = 0.087; // light hole mass
    double mhh = 0.62;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
#endif
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    double mlh = 0.08; // light hole mass
    double mhh = 0.6;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    double mlh = 0.05; // light hole mass
    double mhh = 0.54; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater=="ingap")
  {
    // dnva =  0.5423;
    //mass = 0.665;
    mass = 6.65007285e-1;
  }
  else if (mater == "inp")
  {
    double mlh = 0.074; // light hole mass
    double mhh = 0.5;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else
  {
    Report::UserFatal0() << material << " material not recognized.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::get_DOS_EffectiveMassN
// Purpose       : returns effective mass for electrons for density of states
// Special Notes : See http://ecee.colorado.edu/~bart/book/effmass.htm#dosmass.
//
//   m_e_dos = Mc^(2/3) * (m_l*m_t*m_t)^(1/3)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::get_DOS_EffectiveMassN (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass;

  if (mater == "si")
  {
    double ml = 0.98; // longitudinal mass
    double mt = 0.19; // transverse mass
    double Mc = 6.0; // degeneracy factor (number of equivalent band minimums)
    mass = pow(Mc,2.0/3.0)*pow((ml*mt*mt),1.0/3.0);
    // Note, this should evaluate to 1.08.
  }
  else if (mater == "ge")
  {
    double ml = 1.64; // longitudinal mass
    double mt = 0.082; // transverse mass
    double Mc = 4.0; // degeneracy factor 
    mass = pow(Mc,2.0/3.0)*pow((ml*mt*mt),1.0/3.0);
    // this mass should be around 0.56
  }
  else if (mater == "gaas")
  {
    mass = 0.067; // GaAs is isotropic
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    //dnco=  0.020; = mass^(1.5)
    mass = 0.074;
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    // dnco=  0.0083; = mass^(1.5)
    mass = 0.041;
  }
  else if (mater=="ingap")
  {
    // dnco =  0.02391;
    //mass = 0.0179;
    mass = 8.29952143e-2;
  }
  else if (mater == "inp")
  {
    mass = 0.079;
  }
  else
  {
    Report::UserFatal0() << material << " material not recognized.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::get_DOS_EffectiveMassP
// Purpose       : returns effective mass for holes for density of states
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::get_DOS_EffectiveMassP (const std::string & material)
{
  // Unlike for electrons, Mc is not applied to holes, so this is simply 
  // the same as the non-DOS version.
  return getEffectiveMassP (material);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNc
// Purpose       : density of states, conduction band
// Special Notes :
//
//    Nc = 2*(2*pi*me*kT/h^2)^(3/2)
//       = 2*(2*pi*me_DOS*m0*kT/h^2)^(3/2)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/10/2014
// ----------------------------------------------------------------------------
double MaterialSupport::getNc (const std::string & material, double temp)
{
  double h_planck(CONSTplanck); // Planck's constant  (in J-s)
  double e_mass (CONSTemass);   // e- mass in kg.
  double kb (1.3806488e-23); // boltzmann's constant (J/K)
  double meDOS = get_DOS_EffectiveMassN(material);
  double dnbnd0 = 2.0*M_PI*meDOS*e_mass*kb*temp/(h_planck*h_planck);
  double Nc = 2.0*pow(dnbnd0,1.5)/1.0e6;
  return Nc;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNv
// Purpose       : Density of states, valance band.
// Special Notes :
//
//     Nv = 2*(2*pi*mh*kT/h^2)^(3/2)
//       = 2*(2*pi*mh_DOS*m0*kT/h^2)^(3/2)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/10/2014
// ----------------------------------------------------------------------------
double MaterialSupport::getNv (const std::string & material, double temp)
{
  double h_planck(CONSTplanck); // Planck's constant  (in J-s)
  double e_mass (CONSTemass);   // e- mass in kg.
  double kb (1.3806488e-23); // boltzmann's constant (J/K)
  double mhDOS = get_DOS_EffectiveMassP(material);
  double dnbnd0 = 2.0*M_PI*mhDOS*e_mass*kb*temp/(h_planck*h_planck);
  double Nv = 2.0*pow(dnbnd0,1.5)/1.0e6;
  return Nv;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/12
// ----------------------------------------------------------------------------
double MaterialSupport::getNi (const std::string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;
  double kbq = 8.6173324e-5; // boltzmann's constant  (eV/K)
  double bg = bandgap(mater,temp);
  double Nc = getNc(mater,temp);
  double Nv = getNv(mater,temp);
  ni = sqrt (Nc * Nv) * exp(-bg/(2.0 * kbq * temp));
  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi_old
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getNi_old (const std::string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;

  if (mater == "si")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * pow(6.0,0.5) * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
    // ni = 1.25e10;
  }
  else if (mater == "gaas")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
    //ni = 2.0e6;
  }
  else if (mater == "ge")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * 2.0 * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
   // ni = 2.5e13;
  }
  // for the next several, as they are all III-V materials, I copied the
  // gaas functions.  I *think* this is correct, as I *think* that Mc is
  // going to be 1.0 for all of these.
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater == "inp")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else
  {
    std::string msg = "MaterialSupport::getNi:  ";
    msg += material;
    msg += " material not recognized.\n";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getRelPerm
// Purpose       : returns relative permitivity
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getRelPerm (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double perm;
  if (mater == "si")
  {
    perm = 11.8;
  }
  else if (mater == "sio2")
  {
    perm = 3.9;
  }
  else if (mater == "ge" )
  {
    perm = 16.0;
  }
  else if (mater == "gaas")
  {
    perm = 13.2;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    perm = 12.5;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    perm = 14.0;
  }
  else if (mater == "inp")
  {
    perm = 12.6;
  }
  else
  {
    Report::UserFatal0() << material << " material not recognized.";
  }

  return perm;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRsrh
// Purpose       : Calculates schockley-read-hall recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRsrh
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pn = Ni*Ni;

  double A = (n*p-pn);
  double B = (tp*(n+Ni)+tn*(p+Ni));

  double arg = CONSTMAX_EXP_ARG;
  if (B >= exp(arg)) B = exp(arg);

  return (A/B);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhN
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to electron density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhN
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhN;
  double A1, B1, C1;
  double dAdn;
  double dBdn;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdn = -1.0/(C1*C1) * tp;

  pdRsrhN = dAdn * B1 + dBdn * A1;

  return pdRsrhN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhP
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to hole density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhP
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhP;
  double A1, B1, C1;
  double dAdp;
  double dBdp;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdp = -1.0/(C1*C1) * tn;

  pdRsrhP = dAdp * B1 + dBdp * A1;

  return pdRsrhP;
}


// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRaug
// Purpose       : Calculates Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRaug
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;

  double A = (n*p-pn);
  double C = (Cn*n+Cp*p);

  double arg = CONSTMAX_EXP_ARG;
  if (C >= exp(arg)) C = exp(arg);

  return (A*C);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugN
// Purpose       : Calculates partial derivative w.r.t. electron density
//                 for Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugN
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugN;
  double A1, B1;
  double dAdn;
  double dBdn;

  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdn = Cn;

  pdRaugN = dAdn*B1 + A1*dBdn;

  return pdRaugN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugP
// Purpose       : Calculates partial derivative w.r.t. hole density
//                 for Auger recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 I believe (but am not sure) that the constants Cn and Cp
//                 are material dependent.  That is part of why the
//                 material name is passed in as an argument.  The values
//                 here are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugP
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugP;
  double A1, B1;
  double dAdp;
  double dBdp;

  double Cn = 2.8e-31;
  double Cp = 1.2e-31;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdp = Cp;

  pdRaugP = dAdp*B1 + A1*dBdp;

  return pdRaugP;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::workfunc
// Purpose       : This function returns the workfunction
//                 of various metals
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//----------------------------------------------------------------------------
double MaterialSupport::workfunc(std::string & metal)
{
  double wkfunc=0.0;

  ExtendedString metalName = metal;
  metalName.toLower ();

  if (metalName=="al")
  {
    wkfunc = 4.10;   //aluminum
  }
  else if (metalName=="ppoly")
  {
    wkfunc = 5.25;   // p+-polysilicon
  }
  else if (metalName=="npoly")
  {
    wkfunc = 4.17;   // n+-polysilicon
  }
  else if (metalName=="mo")
  {
    wkfunc = 4.53;  // molybdenum
  }
  else if (metalName=="w")
  {
    wkfunc = 4.63;  // tungsten
  }
  else if (metalName=="modi")
  {
    wkfunc = 4.80;  // molybdenum disilicide
  }
  else if (metalName=="wdi")
  {
    wkfunc = 4.80;  // tungsten disilicide
  }
  else if (metalName=="cu")
  {
    wkfunc = 4.25;   // copper
  }
  else if (metalName=="pt")
  {
    wkfunc = 5.30;   // platinum
  }
  else if (metalName=="au")
  {
    wkfunc = 4.80;   // gold
  }
  else if (metalName=="neutral")
  {
    wkfunc = 0.0;
  }
  else
  {
    Report::UserFatal0() << metalName << " material not recognized.";
  }

  return wkfunc;
}
//----------------------------------------------------------------------------
// Function      : MaterialSupport::affin
// Purpose       : This function returns the electron affinity
//                 of various semiconductor materials
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//---------------------------------------------------------------------------
double MaterialSupport::affin(const std::string & material)
{

  double afty=0.0;

  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")
  {
    afty = 4.17;      // silicon
  }
  else if (materialName=="ge")
  {
    afty = 4.00;     // germanium
  }
  else if (materialName=="gaas")
  {
    afty = 4.07;    // gallium arsenide
  }
  else if (materialName=="sio2")
  {
    afty = 0.97;    // silicon dioxide
  }
  else if (materialName=="nitride")
  {
    afty = 0.97;    // silicon nitride
  }
  else if (materialName=="sapphire")
  {
    afty = 0.97;     // sapphire (also known as aluminum oxide)
  }
  else
  {
    Report::UserError0() << materialName << " material not recognized.";
  }

  return afty;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::bandgap
// Purpose       : This function returns the electronic bandgap
//                 of various semiconductor materials.
//
// Special Notes : Reference for temperature-dependent semiconductor
//                 materials is "The Standard Thermodynamic Function
//                 of the Formation of Electrons and Holes in Ge, Si,
//                 GaAs, and GaP," by C. D. Thurmond, J. Electrochem. Soc.,
//                 vol. 122, p. 1133, 1975.
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/18/03
//---------------------------------------------------------------------------
double MaterialSupport::bandgap(const std::string & material,  double temp)
{
  double gap = 0.0;
  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")   // silicon
  {
    gap = 1.17 - 4.73e-4*pow(temp,2.0)/(temp + 636.0);
  }
  else if (materialName=="ge")  // germanium
  {
    gap = 0.7437 - 4.774e-4*pow(temp,2.0)/(temp + 235);
  }
  else if (materialName=="gaas")  // gallium arsenide
  {
    gap = 1.519 - 5.405e-4*pow(temp,2.0)/(temp + 204);
  }
  else if (materialName=="ingap")  
  {
    gap = 1.86098;
  }
  else if (materialName=="sio2")  // silicon dioxide
  {
    gap = 9.00;
  }
  else if (materialName=="nitride")  // silicon nitride
  {
    gap = 4.7;
  }
  else if (materialName=="sapphire")  // sapphire
  {
    gap = 4.7;
  }
  else if (materialName=="inalas" || materialName=="alinas") // indium aluminum arsenide
  {
    gap = 1.46;
  }
  else if (materialName=="ingaas" || materialName=="gainas") // indium galium arsenide
  {
    gap = 0.75;
  }
  else if (materialName=="inp")  // indium phosphide
  {
    gap = 1.07;
  }
  else
  {
    Report::UserError0() << materialName << " material not recognized.";
  }

  return gap;
}


//----------------------------------------------------------------------------
// Function      : MaterialSupport::Ebgn
// Purpose       : Band-gap narrowing
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::Ebgn (
    const std::string & material,  
    const std::string & bgnModel,
    double dope)
{
  double Ebgn=0.0;
  double n0_bgn=0.0;
  double v0_bgn=0.0;
  double con_bgn=0.0;

  ExtendedString bgnModelName = bgnModel;
  bgnModelName.toLower();

  if (bgnModelName=="slotboom") // parameters taken from Medici manual
  {
    ExtendedString materialName = material;
    materialName.toLower();

    if (materialName=="si")   // silicon
    {
      // Old Slotboom model:
      n0_bgn = 1.0e17; // Nref
      v0_bgn = 9.0e-3; // Eref
      con_bgn = 0.5;
#if 0
      // New Slotboom model:
      double n0_bgn = 1.3e17; // Nref
      double v0_bgn = 6.92e-3; // Eref
      double con_bgn = 0.5;
#endif
    }
    else if (materialName=="ge")  // germanium
    {
      n0_bgn = 1.0e17;
      v0_bgn = 9.0e-3;
      con_bgn = 0.5;
    }
    else if (materialName=="gaas")  // gallium arsenide
    {
      n0_bgn = 1.0e17;
      v0_bgn = 0.0;
      con_bgn = 0.0;
    }
    else if (materialName=="ingap")
    {
      n0_bgn = 1.0e17;
      v0_bgn = 0.0;
      con_bgn = 0.0;
    }
    else if (materialName=="sio2")  // silicon dioxide
    {
      n0_bgn = 1.0e17; // Nref
      v0_bgn = 9.0e-3; // Eref
      con_bgn = 0.5;
    }
#if 0
    else if (materialName=="nitride")  // silicon nitride
    {
      n0_bgn = 1.0;
      v0_bgn = 0.0;
      con_bgn = 0.0;
    }
    else if (materialName=="sapphire")  // sapphire
    {
      n0_bgn = 1.0;
      v0_bgn = 0.0;
      con_bgn = 0.0;
    }
    else if (materialName=="inalas" || materialName=="alinas") // indium aluminum arsenide
    {
      n0_bgn = 1.0e17;
      v0_bgn = 0.0;
      con_bgn = 0.0;
    }
    else if (materialName=="ingaas" || materialName=="gainas") // indium galium arsenide
    {
      Ebgn = 0.12;
    }
    else if (materialName=="inp")  // indium phosphide
    {
      Ebgn = 9.4e-2;
    }
#endif
    else
    {
      Report::UserError0() << materialName << " material not recognized.";
    }

    double tmp1 = log(dope/n0_bgn);
    Ebgn = (v0_bgn) * (tmp1 + sqrt(tmp1*tmp1 + con_bgn));
  }
  else if ("bennet-wilson") // default model for Sentaurus.
  {
    ExtendedString materialName = material;
    materialName.toLower();

    if (materialName=="si") 
    {
      double Eref = 6.84e-3;
      double Nref = 3.162e+18;
      if (dope >= Nref) // Bennet-Wilson model
      {
        Ebgn = Eref*pow(log(dope/Nref),2.0);
      }
    }
    else
    {
      // bennet-wilson not implemented for this material
      Ebgn=0.0;
    }
  }
  else if (bgnModelName=="jain") // consistent with Wampler's XPD code.
  {
    Ebgn = jainEbgn( material,  dope) ;
  }
  else // assume this is the default
  {
    Ebgn = jainEbgn( material,  dope) ;
  }

  return Ebgn;
}


//----------------------------------------------------------------------------
// Function      : MaterialSupport::jainEbgn
//
// Purpose       : Band-gap narrowing, using the Jain model.
//
//                 This corresponds to the BGN2 model in Davinci/Medici
//
// Special Notes : Ntype vs. Ptype is determined by the sign of "dope".
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::jainEbgn (
    const std::string & material,  
    double dope)
{
  double Ebgn=0.0;

  double deltaEcn=0.0;
  double deltaEvn=0.0;
  double deltaEcp=0.0;
  double deltaEvp=0.0;

  double anc_bgn = 0.0;
  double bnc_bgn = 0.0;
  double cnc_bgn = 0.0;
  double anv_bgn = 0.0;
  double bnv_bgn = 0.0;
  double cnv_bgn = 0.0;
  double apc_bgn = 0.0;
  double bpc_bgn = 0.0;
  double cpc_bgn = 0.0;
  double apv_bgn = 0.0;
  double bpv_bgn = 0.0;
  double cpv_bgn = 0.0;

  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")   // silicon
  {  
    anc_bgn = -14.84e-3;
    bnc_bgn = 0.0;
    cnc_bgn = 0.78e-3;

    anv_bgn = 0.0;
    bnv_bgn = 15.08e-3;
    cnv_bgn = 0.74e-3;

    apc_bgn = 0.0;
    bpc_bgn = -16.27e-3;
    cpc_bgn = -0.18e-3;

    apv_bgn = 18.46e-3;
    bpv_bgn = 0.0;
    cpv_bgn = -2.63e-3;
  }
  else if (materialName=="ge")  // germanium
  {
    anc_bgn = -8.67e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -2.02e-3;
    anv_bgn = 0.0;
    bnv_bgn = 8.14e-3;
    cnv_bgn = 2.29e-3;
    apc_bgn = -8.21e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -2.19e-3;
    apv_bgn = 0.0;
    bpv_bgn = 9.18e-3;
    cpv_bgn = 3.58e-3; 
  }
  else if (materialName=="gaas")  // gallium arsenide
  {
    anc_bgn = -16.30e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (materialName=="ingap")
  {
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (materialName=="sio2")  // silicon dioxide (use same as Si)
  {
    anc_bgn = -14.84e-3;
    bnc_bgn = 0.0;
    cnc_bgn = 0.78e-3;

    anv_bgn = 0.0;
    bnv_bgn = 15.08e-3;
    cnv_bgn = 0.74e-3;

    apc_bgn = 0.0;
    bpc_bgn = -16.27e-3;
    cpc_bgn = -0.18e-3;

    apv_bgn = 18.46e-3;
    bpv_bgn = 0.0;
    cpv_bgn = -2.63e-3; 
  }
#if 0
  else if (materialName=="nitride")  // silicon nitride
  {
     // zero
  }
  else if (materialName=="sapphire")  // sapphire
  {
     // zero
  }
#endif
  else if (materialName=="inalas" || materialName=="alinas") // indium aluminum arsenide
  {
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (materialName=="ingaas" || materialName=="gainas") // indium galium arsenide
  {
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (materialName=="inp")  // indium phosphide
  {
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else
  {
    Report::UserError0() << materialName << " material not recognized.";
  }

  double dopeScaled=fabs(dope)/1.0e+18;

  if (dope > 0.0) // n-type
  {
    deltaEcn=
      +anc_bgn*pow(dopeScaled,(1.0/3.0))
      +bnc_bgn*pow(dopeScaled,(0.25))
      +cnc_bgn*pow(dopeScaled,(0.5));

    deltaEvn=
      +anv_bgn*pow(dopeScaled,(1.0/3.0))
      +bnv_bgn*pow(dopeScaled,(0.25))
      +cnv_bgn*pow(dopeScaled,(0.5));

    Ebgn = deltaEcn-deltaEvn;
  }
  else // p-type
  {
    deltaEcp=
      +apc_bgn*pow(dopeScaled,(1.0/3.0))
      +bpc_bgn*pow(dopeScaled,(0.25))
      +cpc_bgn*pow(dopeScaled,(0.5));

    deltaEvp=
      +apv_bgn*pow(dopeScaled,(1.0/3.0))
      +bpv_bgn*pow(dopeScaled,(0.25))
      +cpv_bgn*pow(dopeScaled,(0.5));

    Ebgn = deltaEcp-deltaEvp;
  }


  return Ebgn;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcLt
// Purpose       : This function calculates carrier lifetimes.
//
// Special Notes : holeFlag parameter indicates electrons or holes:
//                     holeFlag = true  -> holes
//                     holeFlag = false -> electrons
//
// This function assumes that conc is an absolute value.
//
// This function comes from this paper:
//
//     "Analysis of High-Efficiency Silicon Solar Cells",
//     IEEE Transactions on Electron Devices, by Harry T.
//     Weaver and R. D. Nasby, vol. ED-28, no. 5, May 1981.
//
// This is a function that probably should have some material dependence,
// but at the moment it doesn't.  I think all the values in this function
// are for Si.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/20/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcLt (bool holeFlag, double conc)
{
  double lt = 0.0;
  double LT0, Nref;

  conc = fabs(conc);

  if (holeFlag)
  {
    LT0   = 3.52e-5;
    Nref  = 7.1e15;
    lt = LT0 / (1.0 + conc / Nref);
  }
  else
  {
    LT0   = 3.95e-4;
    Nref  = 7.1e15;
    lt = LT0 / (1.0 + conc / Nref);
  }
  return lt;
}

} // namespace Device
} // namespace Xyce
