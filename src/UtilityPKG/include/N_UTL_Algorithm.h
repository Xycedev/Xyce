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
// Filename       : $RCSfile: N_UTL_Algorithm.h,v $
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
// Revision Number: $Revision: 1.5.2.1 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_Algorithm_h
#define Xyce_N_UTL_Algorithm_h

#include <string>

#include <N_DEV_fwd.h>
#include <N_DEV_InstanceName.h>

namespace Xyce {
namespace Util {

inline
Device::InstanceName entityNameFromFullParamName(const std::string &full_param_name)
{
  std::string::size_type pos = full_param_name.find(':');

  if (pos == std::string::npos)
    return Device::InstanceName(full_param_name);
  else
    return Device::InstanceName(std::string(full_param_name, 0, pos));
}

inline
std::string paramNameFromFullParamName(const std::string &full_param_name)
{
  std::string::size_type pos = full_param_name.find(':');

  if (pos == std::string::npos)
    return std::string();
  else
    return std::string(full_param_name, pos + 1);
}

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Algorithm_h
