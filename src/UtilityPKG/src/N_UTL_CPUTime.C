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
// Filename       : $RCSfile: N_UTL_CPUTime.C,v $
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
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_CPUTime.h>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#define Xyce_USE_WINDOWS_TIMER
#elif defined(HAVE_SYS_RESOURCE_H)
#include <unistd.h>
#include <sys/resource.h>

#else
#error "Unable to define cpu_time() for an unknown OS."
#endif

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : cpu_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:04:21 2014
//-----------------------------------------------------------------------------
double cpu_time()
{

#if defined(Xyce_USE_WINDOWS_TIMER)
  /* Windows -------------------------------------------------- */
  FILETIME createTime;
  FILETIME exitTime;
  FILETIME kernelTime;
  FILETIME userTime;
  if ( GetProcessTimes( GetCurrentProcess( ),
                        &createTime, &exitTime, &kernelTime, &userTime ) != -1 )
  {
    SYSTEMTIME userSystemTime;
    if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
      return (double)userSystemTime.wHour * 3600.0 +
        (double)userSystemTime.wMinute * 60.0 +
        (double)userSystemTime.wSecond +
        (double)userSystemTime.wMilliseconds / 1000.0;
  }

#elif defined(RUSAGE_SELF)
  {
    struct rusage rusage;
    if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
      return (double) rusage.ru_utime.tv_sec +
        (double) rusage.ru_utime.tv_usec / 1000000.0;
  }
#endif

  return -1.0;		/* Failed. */
}

} // namespace Xyce
