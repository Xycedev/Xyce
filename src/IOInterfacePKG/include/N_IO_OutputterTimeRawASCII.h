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
// Filename       : $RCSfile: N_IO_OutputterTimeRawASCII.h,v $
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
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterTimeRawASCII_h
#define Xyce_N_IO_OutputterTimeRawASCII_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class TimeRawAscii : public TimeInterface
{
public:
  TimeRawAscii(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~TimeRawAscii();

private:
  TimeRawAscii(const TimeRawAscii &);
  TimeRawAscii &operator=(const TimeRawAscii &);

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

  virtual void doStartStep(int step, int max_step)
  {}

  virtual void doSteppingComplete()
  {}

private:
  void timeHeader(Parallel::Machine comm);

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  int                   numPoints_;
  long                  numPointsPos_;
  std::ostream *        os_;
  bool                  printAll_;
  bool                  outputRAWTitleAndDate_;

  Util::Op::OpList      opList_;
};


} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterTimeRawASCII_h
