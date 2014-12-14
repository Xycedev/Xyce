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
// Filename       : $RCSfile: N_DEV_Op.C,v $
//
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2014/08/08 20:02:38 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Op.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceMgrGlobalParameterOp::get
// Purpose       : get the current value of a global param from the device
//                 package
// Special Notes : It is inappropriate for the Device package to be in charge
//                 of global params, but that's where they are right now.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceMgrGlobalParameterOp::get(const DeviceMgrGlobalParameterOp &op, const Util::Op::OpData &op_data)
{
  return op.deviceManager_.getGlobalPar(op.deviceParameterName_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntityParameterOp::get
// Purpose       : get the current value of a device parameter from a device
//                 entity
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceEntityParameterOp::get(const DeviceEntityParameterOp &op, const Util::Op::OpData &op_data)
{
  double result;

  const_cast<DeviceEntity &>(op.deviceEntity_).getParam(op.deviceParameterName_, result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgrParameterOp::get
// Purpose       : get the current value of a device parameter from the device
//                 package
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
DeviceMgrParameterOp::get(const DeviceMgrParameterOp &op, const Util::Op::OpData &op_data)
{
  return op.deviceManager_.getParamNoReduce(op.deviceParameterName_);
}

} // namespace Device
} // namespace Xyce
