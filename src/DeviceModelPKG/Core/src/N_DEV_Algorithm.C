//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Algorithm.C,v $
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/05/19 15:49:14 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Algorithm.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>

#include <N_DEV_ADC.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:31 2014
//-----------------------------------------------------------------------------
///
/// @param device_instance device instance to return name of
///
/// @return name of device instance
///
///
template<>
const std::string &getName(const DeviceInstance *device_instance) 
{
  return device_instance->getName().getEncodedName();
}


//-----------------------------------------------------------------------------
// Function      : getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:31 2014
//-----------------------------------------------------------------------------
///
/// @param device_model device model to return name of
///
/// @return name of device model
///
///
template<>
const std::string &getName(const DeviceModel *device_model) 
{
  return device_model->getName();
}

} // namespace Device
} // namespace Xyce

