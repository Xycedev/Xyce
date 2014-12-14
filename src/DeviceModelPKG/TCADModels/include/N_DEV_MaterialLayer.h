//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, 2013, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MaterialLayer.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/04/30 23:55:34 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MaterialLayer_h
#define Xyce_N_DEV_MaterialLayer_h

#include <string>

#include <N_DEV_CompositeParam.h>
#include <N_DEV_MaterialSupport.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MaterialLayer
// Purpose       :
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 7/14/11
//-----------------------------------------------------------------------------
class MaterialLayer : public CompositeParam
{
  friend class ParametricData<MaterialLayer>;

public:
  static ParametricData<MaterialLayer> &getParametricData();

  MaterialLayer (std::string materialName = std::string("GAAS"));
  virtual ~MaterialLayer ()
  {}

#ifdef Xyce_DEBUG_DEVICE
  friend std::ostream & operator<<(std::ostream & os, const MaterialLayer & ml);
#endif

private:

public:
  std::string name;
  bool nameGiven;
  std::string material;
  bool materialGiven;
  int NX;
  bool NXGiven;
  int LX;
  int begin;  // beginning mesh point
  int end;    // end mesh point +1

  //////////////////////////////////////
  double diel;  // dielectric constant
  bool dielGiven;  // dielectric constant given flag

  double Ec; // conduction band edge
  bool EcGiven; // conduction band edge given flag
  double Ev; // valance band edge
  bool EvGiven ; // valance band edge given flag
  double EcEff; // conduction band edge, including BGN
  double EvEff; // valance band edge, including BGN

  double bg;  // bandgap
  double bgEff; // effective bandgap (including band-gap narrowing, bgn)

  double Cdonor; // n doping concentration
  bool CdonorGiven; // n doping concentration given flag
  double Cacceptor; // p doping concentration
  bool CacceptorGiven; // p doping concentration given flag

  double narco; // band gap narrowing of conduction band
  bool narcoGiven; // band gap narrowing of conduction band given flag
  double narva; // band gap narrowing of valence band
  bool narvaGiven; // band gap narrowing of valence band given flag

  double dnco; // conduction band density of states multiplier = (md*/mo)^3/2
  double dnva;  // valence band density of states multiplier = (md*/mo)^3/2

  double Nc; // conduction band DOS
  double Nv; // valance band DOS

  double emass; // electron DOS effective mass
  bool emassGiven; // electron DOS effective mass given flag
  double hmass; // hole DOS effective mass
  bool hmassGiven; // hole DOS effective mass given flag

  double elmob0; // zero field mobility for electrons (cm2/Vs)
  bool elmob0Given; // zero field mobility for electrons (cm2/Vs) given flag

  double elvsat; // saturation veocity for electrons (cm/s)
  bool elvsatGiven; // saturation veocity for electrons (cm/s) given flag
  double eleo;   // Eo(V/cm) in mobility field dependence

  double homob0; // zero field mobility for holes (cm2/Vs)
  bool homob0Given; // zero field mobility for holes (cm2/Vs) given flag

  double hovsat; // saturation veocity for holes (cm/s)
  bool hovsatGiven; // saturation veocity for holes (cm/s) given flag

  double dir; // direct recombination rate coefficient (cm3/s)

  double augnpp; // Auger recombination rate coefficient for npp (cm3/s)
  double augpnn; // Auger recombination rate coefficient for pnn (cm3/s)

  double srh; // SRH rate coeff (inverse lifetime)
  double srhdet; // energy shift from midgap for SRH

  double Ni; // intrinsic concentration
  bool NiGiven; // intrinsic concentration given flag
  double NiEff; // effective intrinsic concentration, including BGN
  double width;
  bool widthGiven;

  double gradedLayerWidth;
  bool gradedLayerWidthGiven;

  double temperature;

  void processParams ();
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MaterialLayer N_DEV_MaterialLayer;

#endif

