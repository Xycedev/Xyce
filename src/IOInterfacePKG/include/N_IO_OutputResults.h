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
// Filename       : $RCSfile: N_IO_OutputResults.h,v $
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
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputResults_h
#define Xyce_N_IO_OutputResults_h

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace IO {

typedef std::vector<Util::ExpressionData *> ResultVector;

class OutputResults
{
public:
  OutputResults(const std::string &netlist_filename);

  virtual ~OutputResults();

private:
  OutputResults(const OutputResults &);
  OutputResults &operator=(const OutputResults &);

public:
  bool registerPkgOptionsMgr( PkgOptionsMgr &pkgOpt);
  bool addResultParams(const Util::OptionBlock &option_block);

  void output(
    Parallel::Machine                   comm,
    OutputMgr &                         output_manager,
    double                              circuit_time,
    const Analysis::SweepVector &       step_sweep_vector,
    int                                 step_loop_number,
    const N_LAS_Vector &                solution_vector,
    const N_LAS_Vector &                state_vector,
    const N_LAS_Vector &                store_vector);

  void steppingComplete();

private:
  bool                  RESULTinitialized_;
  std::ostream *        os_;                            ///< ostream for .result output
  std::string           netlistFilename_;
  ResultVector          resultVector_;                  ///< Expressions from .RESULT
  bool                  noIndexResult_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputResults_h
