//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_IO_Probe.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/07/09 13:44:48 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <time.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : getTimeDateStamp
// Purpose       : get current date and time and format for .PRINT output
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 06/14/2013
//-----------------------------------------------------------------------------
std::string getTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ", localtime( &now));

  return std::string( timeDate);
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
