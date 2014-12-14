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
// Filename       : $RCSfile: N_DEV_Print.h,v $
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        :
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.4.1 $
//
// Revision Date  : $Date: 2014/08/13 20:36:35 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Print_h
#define Xyce_N_DEV_Print_h

#include <N_DEV_fwd.h>

namespace Xyce {
namespace Device {

std::ostream &print(std::ostream &os, const Device &device);
std::ostream &printDotOpOutput(std::ostream &os, const Device &device);
std::ostream &printParameter(std::ostream &os, const ParameterBase &entity, const std::string &name, const Descriptor &param);
std::ostream &printCompositeParameters(std::ostream &os, const CompositeParam &composite);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Print_h