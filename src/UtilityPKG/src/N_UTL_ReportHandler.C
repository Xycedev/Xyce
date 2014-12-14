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
// Filename       : $RCSfile: N_UTL_ReportHandler.C,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Dave Baur, Raytheon
//
// Creation Date  : 1/7/2014
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/07/31 21:44:39 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <N_UTL_ReportHandler.h>

#include <iostream>
#include <stdexcept>

namespace Xyce {

namespace {

REH s_reportHandler = &default_report_handler;

}

void
report(
  const char *		message,
  unsigned              type)
{
    (*s_reportHandler)(message, type);
}


void
default_report_handler(
  const char *		message,
  unsigned              type)
{
  std::cout << "Message type " << type << ": " << message << std::endl;
}


REH
set_report_handler(
  REH		        reh)
{
  /* %TRACE[ON]% */  /* %TRACE% */
  if (!reh)
    throw std::runtime_error("Cannot set report handler to NULL");

  REH prev_reh = s_reportHandler;
  s_reportHandler = reh;

  return prev_reh;
}

} // namespace Xyce
