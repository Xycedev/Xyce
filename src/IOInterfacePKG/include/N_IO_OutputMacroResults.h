//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_OutputMacroResults.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.1 $
//
// Revision Date  : $Date: 2014/08/28 22:41:47 $
//
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMacroResults_h
#define Xyce_N_IO_OutputMacroResults_h

namespace Xyce {
namespace IO {

// if any macro-simulation level functions were defined, finish them and output results
void outputMacroResults(
  Parallel::Machine             comm,
  const ObjectiveMap &          objective_map,
  const Measure::Manager &      measure_manager,
  FourierMgr &                  fourier_manager,
  std::string                   netlist_filename,
  const StringPairVector &      response_functions,
  std::string                   outputResponseFilename,
  const int                     step_number);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputMacroResults_h
