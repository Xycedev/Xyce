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
// Filename       : $RCSfile: N_UTL_FormatTime.C,v $
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
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_FormatTime.h>

#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : formatTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:16:11 2014
//-----------------------------------------------------------------------------
std::string
formatTime(
  double        time,
  TimeFormat    time_format)
{
  std::stringstream oss;

  if (time < 0.0) {
    time = -time;
    oss << "-";
  }

  if ((time_format & TIMEFORMAT_STYLE_MASK) == TIMEFORMAT_SECONDS) {
    if (time_format & TIMEFORMAT_MILLIS)
      oss << std::fixed << std::setprecision(3) << time;
    else
      oss << std::fixed << std::setprecision(0) << time;
  }
  else if ((time_format & TIMEFORMAT_STYLE_MASK) == TIMEFORMAT_HMS) {
    int int_time = int(time);

    if (time >= 3600.0)
      oss << (int_time)/3600 << ':'
             << std::setw(2) << std::setfill('0') << (int_time/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;

    else if (time >= 60.0)
      oss << ((int) (time)/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;


    else
      oss << ((int) time)%60;

    if (time_format & TIMEFORMAT_MILLIS) {
      int milliseconds = int(std::fmod(time, 1.0)*1000.0 + 0.5);

      oss << '.' << std::setw(3) << std::setfill('0') << milliseconds;
    }
  }
  else
    oss << time;

  return oss.str();
}

} // namespace Xyce
