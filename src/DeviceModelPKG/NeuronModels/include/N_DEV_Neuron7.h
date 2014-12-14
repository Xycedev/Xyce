//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Neuron7.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/05/21 18:25:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron7_h
#define Xyce_N_DEV_Neuron7_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Neuron.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

namespace Xyce {
namespace Device {
namespace Neuron7 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Neuron::Traits>
{
  static const char *name() {return "Neuron";}
  static const char *deviceTypeName() {return "YNEURON level 7";}
  static int numNodes() {return 1;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Neuron device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NeuronInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;
    
public:
  static std::vector< std::vector<int> > jacStamp; // ok for this to be static as all devcies of this type have the same form of jacStamp


  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Miter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  std::map<int,std::string> & getIntNameMap ();
  bool loadDeviceMask ();
  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

public:
  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  Model &       model_;         //< Owning model

  // model level parameters that can be overridden at the instance level
  double memCap;  // membrane capacitance [F]
  double Vt;      // instantaneous threshold potential
  double Vr;      // resting membrane potential
  double Vpeak;   // peak voltage
  double k;       // modeling param []
  double a;       // modeling param
  double b;       // modeling param
  double c;       // modeling param
  double d;       // modeling param
  double uscale;  // u is normally in pico-Amps.  Scalling by 1.0e9 improves accuracy
  double fallRate;// controlls how quickly system resets over the discontinuity

  bool memCapGiven;
  bool VtGiven;
  bool VrGiven;
  bool VpeakGiven;
  bool kGiven;
  bool aGiven;
  bool bGiven;
  bool cGiven;
  bool dGiven;
  bool uscaleGiven;  // u is normally in pico-Amps.  Scalling by 1.0e9 improves accuracy
  bool fallRateGiven;// controlls how quickly system resets over the discontinuity

  // values for loading in the F, Q, dF/dx and dQ/dx
  //
  // equations are:
  // v equation
  //  k * (v-Vr) * (v-Vt) - u + I - C dv/dt = 0
  // u euqtion
  //  a * ( b * (v-Vr) - u ) - du/dt = 0

  double vEquFvalue;
  double vEquQvalue;
  double vEqudFdv;
  double vEqudFdu;
  double vEqudQdv;
  double uEquFvalue;
  double uEquQvalue;
  double uEqudFdv;
  double uEqudFdu;
  double uEqudQdu;

  // these need to be put in the state vector as they
  // make this device state-dependent
  bool resetting;
  double uPeak;

  // offests
  int li_V;
  int li_U;

  // jacobian offsets
  int vEquVOffset;
  int vEquUOffset;
  int uEquVOffset;
  int uEquUOffset;


};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;
    
public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &          MB,
     const FactoryBlock &        factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams ();
  bool processInstanceParams ();

private:
  // model level parameters that can be overridden at the instance level
  double memCap;  // membrane capacitance [F]
  double Vt;      // instantaneous threshold potential
  double Vr;      // resting membrane potential
  double Vpeak;   // peak voltage
  double k;       // modeling param []
  double a;       // modeling param
  double b;       // modeling param
  double c;       // modeling param
  double d;       // modeling param
  double uscale;  // u is normally in pico-Amps.  Scalling by 1.0e9 improves accuracy
  double fallRate;// controlls how quickly system resets over the discontinuity

  bool memCapGiven;
  bool VtGiven;
  bool VrGiven;
  bool VpeakGiven;
  bool kGiven;
  bool aGiven;
  bool bGiven;
  bool cGiven;
  bool dGiven;
  bool uscaleGiven;  // u is normally in pico-Amps.  Scalling by 1.0e9 improves accuracy
  bool fallRateGiven;// controlls how quickly system resets over the discontinuity


public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

private:
  std::vector<Instance*> instanceContainer;
};

void registerDevice();

} // namespace Neuron7
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Neuron7::Instance N_DEV_NeuronInstance7;
typedef Xyce::Device::Neuron7::Model N_DEV_NeuronModel7;

#endif