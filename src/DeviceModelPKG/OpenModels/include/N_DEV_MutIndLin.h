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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MutIndLin.h,v $
//
// Purpose        : Linear Mutual Inductor classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/21/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.75 $
//
// Revision Date  : $Date: 2014/05/22 17:40:29 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MutIndLin_h
#define Xyce_N_DEV_MutIndLin_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : InductorInstancePhysicalData
// Purpose       : This class handles the physical data for an inductor
//                 so that the mutual inductor classes have a clean
//                 contain object.
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
//-----------------------------------------------------------------------------
class InductorInstanceData
{
public:
  InductorInstanceData();
  std::string name;    // name of inductor
  double L;       // User specified inductance
  double IC;      // Initial condition: initial, time-zero inductor current(A)
  bool ICGiven;   // flag if IC was given
  double baseL;   // the baseline inductance before temperature effects
  // local row indicies for external variables
  int li_Pos;     // positive node
  int li_Neg;     // negative node
  int li_Branch;  // branch equation
  // offsets for jacobian enteries
  int APosEquBraVarOffset;
  int ANegEquBraVarOffset;
  int ABraEquPosNodeOffset;
  int ABraEquNegNodeOffset;
  int ABraEquBraVarOffset;
  std::vector< int > inductorCurrentOffsets;
  // offsets for dependent variables
  std::vector<int> ABraEquDepVarOffsets;
  // keep track of which coupling coefficients' dependent vars map to
  // our dependencies.
  std::vector<std::pair<int,int> > depVarPairs;
  // offsets only needed in nonlinear application
  int magOffset;
  int vPosOffset;
  int vNegOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for dFdx entries:
  double * f_PosEquBraVarPtr;
  double * f_NegEquBraVarPtr;
  double * f_BraEquPosNodePtr;
  double * f_BraEquNegNodePtr;
  double * f_BraEquBraVarPtr;
  std::vector< double * > f_inductorCurrentPtrs;
  // offsets for dependent variables
  std::vector<double *> f_BraEquDepVarPtrs;

  // offsets only needed in nonlinear application
  double * f_magPtr;
  double * f_vPosPtr;
  double * f_vNegPtr;

  // Pointers for dQdx entries:
  double * q_PosEquBraVarPtr;
  double * q_NegEquBraVarPtr;
  double * q_BraEquPosNodePtr;
  double * q_BraEquNegNodePtr;
  double * q_BraEquBraVarPtr;
  std::vector< double * > q_inductorCurrentPtrs;
  // offsets for dependent variables
  std::vector<double *> q_BraEquDepVarPtrs;

  // offsets only needed in nonlinear application
  double * q_magPtr;
  double * q_vPosPtr;
  double * q_vNegPtr;
#endif
};

namespace MutIndLin {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Linear Mutual Inductor";}
  static const char *deviceTypeName() {return "K level 1";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "COUP_VAL";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the linear
//                 mutual inductor device.
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Traits;friend class Master;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Iiter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );

  std::map<int,std::string> & getIntNameMap ();
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  const std::vector<std::string> & getDepSolnVars();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);
  void updateInductanceMatrix();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool loadDeviceMask ();
  bool setIC ();

  void varTypes( std::vector<char> & varTypeVec );

  void setupPointers();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void auxDAECalculations ();

  // iterator reference to the inductor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // This container bundles up the physical data for each inductor
  // involved in this mutual inductor
  int numInductors;
  std::vector< InductorInstanceData* > instanceData;

  // These vectors let the new param options load and set inductor data
  // the parser passes all of these to us
  std::vector< std::string > inductorNames;
  std::vector< double > inductorInductances;
  std::vector< std::string > inductorsNode1;
  std::vector< std::string > inductorsNode2;
  // and here's the list of ones we are coupling
  std::vector< std::string > couplingInductor;
  std::vector< double > couplingCoefficient;
  //Pointers to expressions used by coupling coefficients
  std::vector<N_UTL_Expression *> expPtrs;
  // derivatives with respect to expression variables:
  std::vector< std::vector<double> > couplingCoefficientVarDerivs;
  std::vector< std::vector<double> > couplingCoefficientVarDerivsDDT;
  // count of number of variables each coupling coefficient depends on
  std::vector<int> numCoupDepVars;
  std::vector< std::vector< double > > mutualCouplingCoef;
  std::vector< std::vector< double > > mutualCouplingCoefDerivs;
  std::vector<std::pair<int,int> > indexPairs;
  std::vector<int> li_coupling;              // for storage of coupling coeffs in
  // state vector (needed for old DAE)
  std::vector<std::vector<int> > li_couplingVarDerivs; // storing the derivatives
  // of the coupling coefficients
  // for old DAE
  double mutualCup;    // mutual coupling value
  bool mutualCupGiven;

  std::vector< double > inductanceVals;      // the inductances of the inductors
  std::vector< std::vector< double > > LO;        // L' * L (matrix)
  std::vector< double > inductorCurrents;    // currents through inductors (col vec.)
  std::vector< double > dIdt;                // time derivatives of currents.
  // (used by old DAE)
  std::vector< double > LOI;                 // LO * I (col vector).

  double temp;         // temperature of this instance
  bool tempGiven;      // flag if temp was given

  // scaling crontrol for new equations
  double scalingRHS;   // scaling for loading DAE or RHS components

  // this can't be static as each instance may have a different
  // number of inductors in it so they'll each have a different
  // size jacStamp
  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Traits;friend class Master;

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


public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

private:
  std::vector<Instance*> instanceContainer;

private:

  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double tnom;         // reference temperature
  // flags indicating if temperature parameters were given
  bool tc1Given;
  bool tc2Given;
  bool tnomGiven;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // new DAE stuff:
  // new DAE load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * storeLeadF, double * storeLeadQ);
  virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
};

void registerDevice();

} // namespace MutIndLin
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MutIndLin::Instance N_DEV_MutIndLinInstance;
typedef Xyce::Device::MutIndLin::Model N_DEV_MutIndLinModel;
typedef Xyce::Device::MutIndLin::Master N_DEV_MutIndLinMaster;

#endif

