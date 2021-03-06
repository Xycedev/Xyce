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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Neuron9.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Christy Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.27 $
//
// Revision Date  : $Date: 2014/05/13 14:50:39 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Neuron9.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Neuron.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {


namespace Neuron9 {


void Traits::loadInstanceParameters(ParametricData<Neuron9::Instance> &p)
{
}

void Traits::loadModelParameters(ParametricData<Neuron9::Model> &p)
{
  // Set up map for double precision variables:
  p.addPar ("CMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::cMem,
          &Neuron9::Model::cMemGiven,
          U_FARAD, CAT_NONE, "Membrane capacitance");

  p.addPar ("ELEAK", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::eLeak,
          &Neuron9::Model::eLeakGiven,
          U_VOLT, CAT_NONE, "Leak current reversal potential");

  p.addPar ("GMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::gMem,
          &Neuron9::Model::gMemGiven,
          U_OHMM1, CAT_NONE, "Membrane conductance");

  p.addPar ("EK", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::eK,
          &Neuron9::Model::eKGiven,
          U_VOLT, CAT_NONE, "Potassium reversal potential");

  p.addPar ("GK", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::gK,
          &Neuron9::Model::gKGiven,
          U_OHMM1, CAT_NONE, "Potassium base conductance");

  p.addPar ("ENA", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::eNa,
          &Neuron9::Model::eNaGiven,
          U_VOLT, CAT_NONE, "Sodium reversal potential");

  p.addPar ("GNA", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::gNa,
          &Neuron9::Model::gNaGiven,
          U_OHMM1, CAT_NONE, "Sodium base conductance");

  p.addPar ("VREST", 0.0, false, ParameterType::NO_DEP,
          &Neuron9::Model::vRest,
          &Neuron9::Model::vRestGiven,
          U_VOLT, CAT_NONE, "Resting potential");
}



//
// static class member inits
//
std::vector< std::vector<int> > Instance::jacStamp;

const double Instance::VT = -63.0;	// mV

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)

  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    kcl1Fvalue(0.0),
    kcl1Qvalue(0.0),
    kcl2Fvalue(0.0),
    kcl2Qvalue(0.0),
    nEquFvalue(0.0),
    nEquQvalue(0.0),
    mEquFvalue(0.0),
    mEquQvalue(0.0),
    hEquFvalue(0.0),
    hEquQvalue(0.0),
    dkcl1F_dV1(0.0),
    dkcl1F_dV2(0.0),
    dkcl1F_dn(0.0),
    dkcl1F_dm(0.0),
    dkcl1F_dh(0.0),
    dkcl1Q_dV1(0.0),
    dkcl1Q_dV2(0.0),
    dkcl2F_dV1(0.0),
    dkcl2F_dV2(0.0),
    dkcl2F_dn(0.0),
    dkcl2F_dm(0.0),
    dkcl2F_dh(0.0),
    dkcl2Q_dV1(0.0),
    dkcl2Q_dV2(0.0),
    dnF_dV1(0.0),
    dnF_dV2(0.0),
    dnF_dn(0.0),
    dnQ_dn(0.0),
    dmF_dV1(0.0),
    dmF_dV2(0.0),
    dmF_dm(0.0),
    dmQ_dm(0.0),
    dhF_dV1(0.0),
    dhF_dV2(0.0),
    dhF_dh(0.0),
    dhQ_dh(0.0)
{
  numExtVars   = 2;
  numIntVars   = 3;
  numStateVars = 2;

  devConMap.resize(2);
  devConMap[0] = 1;
  devConMap[1] = 1;

  // set up jacStamp
  if( jacStamp.empty() )
  {
    // V1, V2, n, m, h
    // all values depend on Vm = V1-V2
    // V1 and V2 depend on n, m, and h
    // n, m, and h each depend on themselves as well as V1 and V2
    jacStamp.resize(5);
    jacStamp[0].resize(5);	// V1
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[0][3] = 3;
    jacStamp[0][4] = 4;
    jacStamp[1].resize(5);	// V2
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[1][3] = 3;
    jacStamp[1][4] = 4;
    jacStamp[2].resize(3);	// n
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
    jacStamp[2][2] = 2;
    jacStamp[3].resize(3);	// m
    jacStamp[3][0] = 0;
    jacStamp[3][1] = 1;
    jacStamp[3][2] = 3;
    jacStamp[4].resize(3);	// h
    jacStamp[4][0] = 0;
    jacStamp[4][1] = 1;
    jacStamp[4][2] = 4;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  // there are no params for this device as yet, so calling setParams causes
  // the code to exit.
  // setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                            const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Pos  = extLIDVec[0];
  li_Neg  = extLIDVec[1];
  li_nPro = intLIDVec[0];
  li_mPro = intLIDVec[1];
  li_hPro = intLIDVec[2];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl
         << "  li_Neg = " << li_Neg << std::endl
         << "  li_nPro = " << li_nPro << std::endl
         << "  li_mPro = " << li_mPro << std::endl
         << "  li_hPro = " << li_hPro << std::endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    intNameMap[ li_nPro ] = spiceInternalName(getName(), "N");
    intNameMap[ li_mPro ] = spiceInternalName(getName(), "M");
    intNameMap[ li_hPro ] = spiceInternalName(getName(), "H");
  }
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  li_KCurrentState = staLIDVec[0];
  li_NaCurrentState = staLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDeviceMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  bool returnVal=false;
  N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;

//   Xyce::dout() << "Masking n, m and h" << std::endl;
//   (*maskVectorPtr)[li_nPro] = 0.0;
//   (*maskVectorPtr)[li_mPro] = 0.0;
//   (*maskVectorPtr)[li_hPro] = 0.0;
//   returnVal = true;

  return (returnVal);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  APosEquNNodeOffset   = jacLIDVec[0][2];
  APosEquMNodeOffset   = jacLIDVec[0][3];
  APosEquHNodeOffset   = jacLIDVec[0][4];

  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquNNodeOffset   = jacLIDVec[1][2];
  ANegEquMNodeOffset   = jacLIDVec[1][3];
  ANegEquHNodeOffset   = jacLIDVec[1][4];

  ANEquPosNodeOffset   = jacLIDVec[2][0];
  ANEquNegNodeOffset   = jacLIDVec[2][1];
  ANEquNNodeOffset     = jacLIDVec[2][2];

  AMEquPosNodeOffset   = jacLIDVec[3][0];
  AMEquNegNodeOffset   = jacLIDVec[3][1];
  AMEquMNodeOffset     = jacLIDVec[3][2];

  AHEquPosNodeOffset   = jacLIDVec[4][0];
  AHEquNegNodeOffset   = jacLIDVec[4][1];
  AHEquHNodeOffset     = jacLIDVec[4][2];
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  //Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl;

  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  // use suffix "now" to to clarify that this for the latest solution
  double v1Now = (*solVectorPtr)[li_Pos];
  double v2Now = (*solVectorPtr)[li_Neg];
  double nNow  = (*solVectorPtr)[li_nPro];
  double mNow  = (*solVectorPtr)[li_mPro];
  double hNow  = (*solVectorPtr)[li_hPro];

  // get function and derivative values
  // independent variables
  // use scoping to avoid lots of similar variable names

  // kcl F
  {
    const int numDeriv = 5;
    Sacado::Fad::SFad<double,5> v1Var( numDeriv, 0, v1Now );
    Sacado::Fad::SFad<double,5> v2Var( numDeriv, 1, v2Now );
    Sacado::Fad::SFad<double,5> nVar( numDeriv, 2, nNow );
    Sacado::Fad::SFad<double,5> mVar( numDeriv, 3, mNow );
    Sacado::Fad::SFad<double,5> hVar( numDeriv, 4, hNow );
    // parameters
    Sacado::Fad::SFad<double,5> gMemVar( model_.gMem );
    Sacado::Fad::SFad<double,5> eLeakVar( model_.eLeak );
    Sacado::Fad::SFad<double,5> gKVar( model_.gK );
    Sacado::Fad::SFad<double,5> eKVar( model_.eK );
    Sacado::Fad::SFad<double,5> gNaVar( model_.gNa );
    Sacado::Fad::SFad<double,5> eNaVar( model_.eNa );

    // compute the vaud and derivative terms for KCL 1 F
    Sacado::Fad::SFad<double,5> resultFad;
    resultFad = kcl1EquF( v1Var, v2Var, nVar, mVar, hVar, gMemVar, eLeakVar, gKVar, eKVar, gNaVar, eNaVar );
    kcl1Fvalue = resultFad.val();
    dkcl1F_dV1 = resultFad.dx(0);
    dkcl1F_dV2 = resultFad.dx(1);
    dkcl1F_dn  = resultFad.dx(2);
    dkcl1F_dm  = resultFad.dx(3);
    dkcl1F_dh  = resultFad.dx(4);

    // compute the vaud and derivative terms for KCL 2 F
    resultFad = kcl2EquF( v1Var, v2Var, nVar, mVar, hVar, gMemVar, eLeakVar, gKVar, eKVar, gNaVar, eNaVar );
    kcl2Fvalue = resultFad.val();
    dkcl2F_dV1 = resultFad.dx(0);
    dkcl2F_dV2 = resultFad.dx(1);
    dkcl2F_dn  = resultFad.dx(2);
    dkcl2F_dm  = resultFad.dx(3);
    dkcl2F_dh  = resultFad.dx(4);
  }

  // kcl Q
  {
    const int numDeriv = 2;
    Sacado::Fad::SFad<double,2> v1Var( numDeriv, 0, v1Now );
    Sacado::Fad::SFad<double,2> v2Var( numDeriv, 1, v2Now );

    // parameters
    Sacado::Fad::SFad<double,2> cMemVar( model_.cMem );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = kcl1EquQ( v1Var, v2Var, cMemVar );
    kcl1Qvalue = resultFad.val();
    dkcl1Q_dV1 = resultFad.dx(0);
    dkcl1Q_dV2 = resultFad.dx(1);

    resultFad = kcl2EquQ( v1Var, v2Var, cMemVar );
    kcl2Qvalue = resultFad.val();
    dkcl2Q_dV1 = resultFad.dx(0);
    dkcl2Q_dV2 = resultFad.dx(1);
  }

  // n - equation
  {
    const int numDeriv = 3;
    Sacado::Fad::SFad<double,3> v1Var( numDeriv, 0, v1Now );
    Sacado::Fad::SFad<double,3> v2Var( numDeriv, 1, v2Now );
    Sacado::Fad::SFad<double,3> nVar( numDeriv, 2, nNow );
    // parameter
    Sacado::Fad::SFad<double,3> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,3> resultFad = nEquF( v1Var, v2Var, nVar, vRestVar);
    nEquFvalue = resultFad.val();
    dnF_dV1 = resultFad.dx(0);
    dnF_dV2 = resultFad.dx(1);
    dnF_dn  = resultFad.dx(2);
  }
  {
    const int numDeriv = 1;
    Sacado::Fad::SFad<double,1> nVar( numDeriv, 0, nNow );

    Sacado::Fad::SFad<double,1> resultFad = nEquQ( nVar );
    nEquQvalue = resultFad.val();
    dnQ_dn     = resultFad.dx(0);
  }

  // m - equation
  {
    const int numDeriv = 3;
    Sacado::Fad::SFad<double,3> v1Var( numDeriv, 0, v1Now );
    Sacado::Fad::SFad<double,3> v2Var( numDeriv, 1, v2Now );
    Sacado::Fad::SFad<double,3> mVar( numDeriv, 2, mNow );
    // parameter
    Sacado::Fad::SFad<double,3> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,3> resultFad = mEquF( v1Var, v2Var, mVar, vRestVar );
    mEquFvalue = resultFad.val();
    dmF_dV1 = resultFad.dx(0);
    dmF_dV2 = resultFad.dx(1);
    dmF_dm  = resultFad.dx(2);
  }
  {
    const int numDeriv = 1;
    Sacado::Fad::SFad<double,1> mVar( numDeriv, 0, mNow );

    Sacado::Fad::SFad<double,1> resultFad = mEquQ( mVar );
    mEquQvalue = resultFad.val();
    dmQ_dm     = resultFad.dx(0);
  }

  // h - equation
  {
    const int numDeriv = 3;
    Sacado::Fad::SFad<double,3> v1Var( numDeriv, 0, v1Now );
    Sacado::Fad::SFad<double,3> v2Var( numDeriv, 1, v2Now );
    Sacado::Fad::SFad<double,3> hVar( numDeriv, 2, hNow );
    // parameter
    Sacado::Fad::SFad<double,3> vRestVar( model_.vRest );

    Sacado::Fad::SFad<double,3> resultFad = hEquF( v1Var, v2Var, hVar, vRestVar );
    hEquFvalue = resultFad.val();
    dhF_dV1 = resultFad.dx(0);
    dhF_dV2 = resultFad.dx(1);
    dhF_dh  = resultFad.dx(2);
  }
  {
    const int numDeriv = 1;
    Sacado::Fad::SFad<double,1> hVar( numDeriv, 0, hNow );

    Sacado::Fad::SFad<double,1> resultFad = hEquQ( hVar );
    hEquQvalue = resultFad.val();
    dhQ_dh     = resultFad.dx(0);
  }

  // calculate potassium current
  potassiumCurrent = model_.gK * pow(nNow, 4.0) * (v1Now - v2Now - model_.eK);

  // calculate sodium current
  sodiumCurrent = model_.gNa * pow(mNow, 3.0) * hNow * (v1Now - v2Now - model_.eNa);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    double vRest = model_.vRest;
    Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl
              << "vRest(input) = " << vRest << std::endl
              << "v1 = " << v1Now << std::endl
              << "v2 = " << v2Now << std::endl
              << "nNow = " << nNow << std::endl
              << "mNow = " << mNow << std::endl
              << "hNow = " << hNow << std::endl
              << "kcl1Fvalue = " << kcl1Fvalue << std::endl
              << "dkcl1F_dV1 = " << dkcl1F_dV1 << std::endl
              << "dkcl1F_dV2 = " << dkcl1F_dV2 << std::endl
              << "dkcl1F_dn = " << dkcl1F_dn << std::endl
              << "dkcl1F_dm = " << dkcl1F_dm << std::endl
              << "dkcl1F_dh = " << dkcl1F_dh << std::endl
              << "kcl2Fvalue = " << kcl2Fvalue << std::endl
              << "dkcl2F_dV1 = " << dkcl2F_dV1 << std::endl
              << "dkcl2F_dV2 = " << dkcl2F_dV2 << std::endl
              << "dkcl2F_dn = " << dkcl2F_dn << std::endl
              << "dkcl2F_dm = " << dkcl2F_dm << std::endl
              << "dkcl2F_dh = " << dkcl2F_dh << std::endl
              << "alphaN = " << alphaN<double>( v1Now, v2Now, vRest) << std::endl
              << "betaN = " << betaN<double>( v1Now, v2Now, vRest) << std::endl
              << "nEquFvalue = " << nEquFvalue << std::endl
              << "dnF_dV1 = " << dnF_dV1 << std::endl
              << "dnF_dV2 = " << dnF_dV2 << std::endl
              << "dnF_dn = " << dnF_dn << std::endl
              << "nEquQvalue = " << nEquQvalue << std::endl
              << "dnQ_dn = " << dnQ_dn << std::endl
              << "alphaM = " << alphaM<double>( v1Now, v2Now, vRest) << std::endl
              << "betaM = " << betaM<double>( v1Now, v2Now, vRest) << std::endl
              << "mEquFvalue = " << mEquFvalue << std::endl
              << "dmF_dV1 = " << dmF_dV1 << std::endl
              << "dmF_dV2 = " << dmF_dV2 << std::endl
              << "dmF_dm = " << dmF_dm << std::endl
              << "mEquQvalue = " << mEquQvalue << std::endl
              << "dmQ_dm = " << dmQ_dm << std::endl
              << "alphaH = " << alphaH<double>( v1Now, v2Now, vRest) << std::endl
              << "betaH = " << betaH<double>( v1Now, v2Now, vRest) << std::endl
              << "hEquFvalue = " << hEquFvalue << std::endl
              << "dhF_dV1 = " << dhF_dV1 << std::endl
              << "dhF_dV2 = " << dhF_dV2 << std::endl
              << "dhF_dh = " << dhF_dh << std::endl
              << "hEquQvalue = " << hEquQvalue << std::endl
              << "dhQ_dh = " << dhQ_dh << std::endl
              << std::endl;
  }
#endif

  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars();

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  (*staVectorPtr)[li_KCurrentState]  = potassiumCurrent;
  (*staVectorPtr)[li_NaCurrentState] = sodiumCurrent;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  N_LAS_Vector * daeQVecPtr = extData.daeQVectorPtr;

  (*daeQVecPtr)[li_Pos]  += kcl1Qvalue;
  (*daeQVecPtr)[li_Neg]  += kcl2Qvalue;
  (*daeQVecPtr)[li_nPro] += nEquQvalue;
  (*daeQVecPtr)[li_mPro] += mEquQvalue;
  (*daeQVecPtr)[li_hPro] += hEquQvalue;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron instance.
//
// Special Notes :
//
//
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_Pos]  += kcl1Fvalue;
  (*daeFVecPtr)[li_Neg]  += kcl2Fvalue;
  (*daeFVecPtr)[li_nPro] += nEquFvalue;
  (*daeFVecPtr)[li_mPro] += mEquFvalue;
  (*daeFVecPtr)[li_hPro] += hEquFvalue;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron instance.
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  (*dQdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1Q_dV1;
  (*dQdxMatPtr)[li_Pos][APosEquNegNodeOffset] += dkcl1Q_dV2;

  (*dQdxMatPtr)[li_Neg][ANegEquPosNodeOffset] += dkcl2Q_dV1;
  (*dQdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2Q_dV2;

  (*dQdxMatPtr)[li_nPro][ANEquNNodeOffset] += dnQ_dn;

  (*dQdxMatPtr)[li_mPro][AMEquMNodeOffset] += dmQ_dm;

  (*dQdxMatPtr)[li_hPro][AHEquHNodeOffset] += dhQ_dh;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1F_dV1;
  (*dFdxMatPtr)[li_Pos][APosEquNegNodeOffset] += dkcl1F_dV2;
  (*dFdxMatPtr)[li_Pos][APosEquNNodeOffset]   += dkcl1F_dn;
  (*dFdxMatPtr)[li_Pos][APosEquMNodeOffset]   += dkcl1F_dm;
  (*dFdxMatPtr)[li_Pos][APosEquHNodeOffset]   += dkcl1F_dh;

  (*dFdxMatPtr)[li_Neg][ANegEquPosNodeOffset] += dkcl2F_dV1;
  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2F_dV2;
  (*dFdxMatPtr)[li_Neg][ANegEquNNodeOffset]   += dkcl2F_dn;
  (*dFdxMatPtr)[li_Neg][ANegEquMNodeOffset]   += dkcl2F_dm;
  (*dFdxMatPtr)[li_Neg][ANegEquHNodeOffset]   += dkcl2F_dh;

  (*dFdxMatPtr)[li_nPro][ANEquPosNodeOffset] += dnF_dV1;
  (*dFdxMatPtr)[li_nPro][ANEquNegNodeOffset] += dnF_dV2;
  (*dFdxMatPtr)[li_nPro][ANEquNNodeOffset]   += dnF_dn;

  (*dFdxMatPtr)[li_mPro][AMEquPosNodeOffset] += dmF_dV1;
  (*dFdxMatPtr)[li_mPro][AMEquNegNodeOffset] += dmF_dV2;
  (*dFdxMatPtr)[li_mPro][AMEquMNodeOffset]   += dmF_dm;

  (*dFdxMatPtr)[li_hPro][AHEquPosNodeOffset] += dhF_dV1;
  (*dFdxMatPtr)[li_hPro][AHEquNegNodeOffset] += dhF_dV2;
  (*dFdxMatPtr)[li_hPro][AHEquHNodeOffset]   += dhF_dh;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

// additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Neuron instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << std::endl;
  }

  os << std::endl;
  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}



//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    // just call the instance level update state
    (*it)->updatePrimaryState();
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("neuron", 9)
    .registerModelType("neuron", 9);
}

} // namespace Neuron9
} // namespace Device
} // namespace Xyce
