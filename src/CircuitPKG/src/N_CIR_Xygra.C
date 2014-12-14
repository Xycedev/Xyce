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
// Filename       : $RCSfile: N_CIR_Xygra.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 8/21/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/04/21 16:58:04 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <N_CIR_Xygra.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_Xygra.h>

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getXygraInstance_
// Purpose       : Returns the pointer to a named Xygra device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Xyce::Device::Xygra::Instance *
N_CIR_Xygra::getXygraInstance_(const std::string & deviceName)
{
  // See if we've looked this up before.
  if (xygraDeviceMap_.empty())
  {
    Xyce::Device::Device *device = devIntPtr_->getDevice(Xyce::Device::Xygra::Traits::modelGroup());
    if (device)
      Xyce::Device::mapDeviceInstances(*device, xygraDeviceMap_);
  }

  std::map<std::string, Xyce::Device::Xygra::Instance *>::iterator mapIter = xygraDeviceMap_.find(deviceName);
  if (mapIter == xygraDeviceMap_.end())
    return 0;

  return (*mapIter).second;
}


//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraGetNumNodes
// Purpose       : Returns the number of nodes that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int N_CIR_Xygra::xygraGetNumNodes(const std::string & deviceName)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->getNumNodes();
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraGetNumWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int N_CIR_Xygra::xygraGetNumWindings(const std::string & deviceName)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->getNumWindings();
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraGetCoilWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has in each coil, given the name of that instance in the
//                 netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/15/08
//-----------------------------------------------------------------------------
void N_CIR_Xygra::xygraGetCoilWindings(const std::string & deviceName,
                                           std::vector<int> & cW)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  xygraInstancePtr->getCoilWindings(cW);
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraGetCoilNames
// Purpose       : Returns the names of each coil in a given Xygra instance
//                 given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/29/08
//-----------------------------------------------------------------------------
void N_CIR_Xygra::xygraGetCoilNames(const std::string & deviceName,
                                           std::vector<std::string> & cN)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  xygraInstancePtr->getCoilNames(cN);
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraSetConductances
// Purpose       : Sets th conductance matrix on the specified Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetConductances(const std::string & deviceName,
                                           const std::vector<std::vector<double> > & cM)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->setConductances(cM);
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraSetK
// Purpose       : Sets the K matrix on the specified Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetK(const std::string & deviceName,
                                const std::vector<std::vector<double> > & kM,
                                const double t)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->setK(kM,t);
}

//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraSetSources
// Purpose       : Set the S vector on named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetSources(const std::string & deviceName,
                                      const std::vector<double> & sV,
                                      const double t)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->setSources(sV,t);
}


//-----------------------------------------------------------------------------
// Function      : N_CIR_Xygra::xygraGetVoltages
// Purpose       : Retrieve the voltages on nodes for named Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraGetVoltages(const std::string & deviceName,
                                       std::vector<double> & vN)
{
  Xyce::Device::Xygra::Instance * xygraInstancePtr = getXygraInstance_(deviceName);
  return xygraInstancePtr->getVoltages(vN);
}

