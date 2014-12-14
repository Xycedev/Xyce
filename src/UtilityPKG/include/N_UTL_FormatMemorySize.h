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
// Filename       : $RCSfile: N_UTL_FormatMemorySize.h,v $
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
// Revision Number: $Revision: 1.1.2.1 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_FormatMemorySize_h
#define Xyce_N_UTL_FormatMemorySize_h

#include <string>

namespace Xyce {

typedef size_t MemorySize;

//-----------------------------------------------------------------------------
// Function      : formatMemorySize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:06:44 2014
//-----------------------------------------------------------------------------
///
/// Returns a string of the memory size.
///
/// @param size memory to stringize
///
/// @return string of memory size
///
std::string formatMemorySize(double size);

//-----------------------------------------------------------------------------
// Function      : formatMemorySize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:06:44 2014
//-----------------------------------------------------------------------------
///
/// Returns a string of the memory size.
///
/// @param size memory to stringize
///
/// @return string of memory size
///
std::string formatMemorySize(MemorySize size);

} // namespace Xyce

#endif // Xyce_N_UTL_FormatMemorySize_h
