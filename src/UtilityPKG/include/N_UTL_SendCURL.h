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
// Filename       : $RCSfile: N_UTL_SendCURL.h,v $
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
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2014/08/26 22:04:40 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_SendCURL_h
#define Xyce_N_UTL_SendCURL_h

#include <string>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : sendTrackingData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Aug 26 15:51:21 2014
//-----------------------------------------------------------------------------
///
/// Send tracking data to the URL specified through the specified proxy.
/// If url if null or empty, no tracking data is sent.
///
/// @param url          pointer to URL to send data; null or empty to
///                     disable sending of data
/// @param proxy        pointer to URL of proxy; null or empty for no proxy
/// @param post_message message string
///
/// @return true if message sent, or if url null or empty
///
bool
sendTrackingData(
  const char *          url,
  const char *          proxy,
  const std::string &   post_message);

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_SendCURL_h
