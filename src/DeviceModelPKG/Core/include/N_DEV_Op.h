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
// Filename       : $RCSfile: N_DEV_Op.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/08/08 20:02:37 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Op_h
#define Xyce_N_DEV_Op_h

#include <iterator>

#include <N_DEV_fwd.h>
#include <N_UTL_Op.h>

namespace Xyce {
namespace Device {

class DeviceMgrGlobalParameterOp : public Util::Op::Op<DeviceMgrGlobalParameterOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  DeviceMgrGlobalParameterOp(const std::string &name, const DeviceMgr &device_manager, const std::string &device_parameter_name)
    : Base(name),
      deviceManager_(device_manager),
      deviceParameterName_(device_parameter_name)
  {}

  virtual ~DeviceMgrGlobalParameterOp()
  {}

  static complex get(const DeviceMgrGlobalParameterOp &op, const Util::Op::OpData &op_data);

  const std::string             deviceParameterName_;
  const DeviceMgr &             deviceManager_;
};

class DeviceEntityParameterOp : public Util::Op::Op<DeviceEntityParameterOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  DeviceEntityParameterOp(const std::string &name, const DeviceEntity &device_entity, const std::string &device_parameter_name)
    : Base(name),
      deviceEntity_(device_entity),
      deviceParameterName_(device_parameter_name)
  {}
    

  virtual ~DeviceEntityParameterOp()
  {}

  static complex get(const DeviceEntityParameterOp &op, const Util::Op::OpData &op_data);

  const std::string             deviceParameterName_;
  const DeviceEntity &          deviceEntity_;
};

class DeviceMgrParameterOp : public Util::Op::Op<DeviceMgrParameterOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  DeviceMgrParameterOp(const std::string &name, const DeviceMgr &device_manager, const std::string &device_parameter_name)
    : Base(name),
      deviceManager_(device_manager),
      deviceParameterName_(device_parameter_name)
  {}


  virtual ~DeviceMgrParameterOp()
  {}

  static complex get(const DeviceMgrParameterOp &op, const Util::Op::OpData &op_data);

  const std::string     deviceParameterName_;
  const DeviceMgr &     deviceManager_;
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Op_h
