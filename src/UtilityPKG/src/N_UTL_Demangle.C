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
// Filename       : $RCSfile: N_UTL_Demangle.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 3/20/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/03/19 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#include <N_UTL_Demangle.h>
#include <stdlib.h>

#if __GNUC__ == 3 || __GNUC__ == 4
#include <cxxabi.h>
#endif

// #if defined __xlC__
// #include <demangle.h>
// #endif

namespace Xyce {

#ifdef Xyce__USE_PLATFORM_DEMANGLER

#if defined(__GNUC__)

#if (__GNUC__ == 3)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;
  int		status = 0;

  char *demangled_symbol = abi::__cxa_demangle(symbol, 0, 0, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}

#elif (__GNUC__ == 4)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;

  int		status;

  char *demangled_symbol = __cxxabiv1::__cxa_demangle(symbol, 0, 0, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}

#endif // (__GNUC__ == 3)

#elif defined __xlC__
std::string
demangle(
  const char *	symbol)
{
  return symbol;
// #ifdef PURIFY_BUILD
//   return symbol;
// #else
//   char *rest;

//   Name *name = Demangle(symbol, rest) ;

//   std::string s(name ? name->Text() : symbol);

//   delete name;

//   return s;
// #endif
}

#endif // defined __GNUC__
#else
std::string demangle(const char *symbol) 
{
  return symbol;
}
#endif // Xyce__USE_PLATFORM_DEMANGLER

} // namespace Xyce
