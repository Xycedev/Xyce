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
// Filename       : $RCSfile: N_UTL_WallTime.C,v $
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
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_WallTime.h>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#define Xyce_USE_WINDOWS_TIMER
#else
#include <sys/time.h>
#endif

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : wall_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:34:46 2014
//-----------------------------------------------------------------------------
double
wall_time()
{
#ifdef Xyce_USE_WINDOWS_TIMER
  LARGE_INTEGER time,freq;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&time);
  return static_cast<double>(time.QuadPart)/static_cast<double>(freq.QuadPart);
#else
  timeval tp;
  struct timezone tz;
  ::gettimeofday(&tp, &tz);

  double seconds = tp.tv_sec;
  double milliseconds = tp.tv_usec*1.0e-6;

  return seconds + milliseconds;
#endif
}

//-----------------------------------------------------------------------------
// Function      : wall_dtime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:34:54 2014
//-----------------------------------------------------------------------------
double
wall_dtime(double &t)
{
  const double tnew = wall_time();

  const double dt = tnew - t;

  t = tnew ;

  return dt ;
}

} // namespace Xyce
