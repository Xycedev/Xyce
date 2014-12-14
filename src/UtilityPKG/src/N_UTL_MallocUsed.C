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
// Filename       : $RCSfile: N_UTL_MallocUsed.C,v $
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
// Revision Number: $Revision: 1.1.2.4 $
//
// Revision Date  : $Date: 2014/09/03 16:40:08 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <cstddef>

#if defined(HAVE_MALLOC_H) && defined(HAVE_MALLINFO)
#include <malloc.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Function      : malloc_used
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:55:07 2014
//-----------------------------------------------------------------------------
size_t
malloc_used()
{
  size_t heap_size = 0;
  size_t largest_free = 0;

#if defined(HAVE_MALLOC_H) && defined(HAVE_MALLINFO)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = (size_t) minfo.uordblks + (size_t) minfo.hblkhd;
  largest_free = (size_t) minfo.fordblks;
#endif

  return heap_size;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
