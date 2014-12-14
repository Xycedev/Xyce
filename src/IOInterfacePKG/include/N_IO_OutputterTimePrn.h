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
// Filename       : $RCSfile: N_IO_OutputterTimePrn.h,v $
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
// Revision Number: $Revision: 1.7.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterTimePrn_h
#define Xyce_N_IO_OutputterTimePrn_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : TimePrn
// Purpose       : Outputter class for transient runs and "PRN" (Standard)
//                 output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
class TimePrn : public TimeInterface
{
public:
  static Outputter::Interface *factory(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);
  static void registerOutputter(OutputMgr &output_manager);

  TimePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~TimePrn();

private:
  TimePrn(const TimePrn &);
  TimePrn &operator=(const TimePrn &);

public:
  virtual void doSetAnalysisMode(Analysis::Analysis_Mode analysis_mode) {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputTime(
    Parallel::Machine           comm,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector) ;

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doSteppingComplete();

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  std::ostream *        os_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      opList_;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterTimePrn_h
