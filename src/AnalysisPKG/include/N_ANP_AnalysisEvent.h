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
// Filename       : $RCSfile: N_ANP_AnalysisEvent.h,v $
//
// Purpose        : Step analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2014/08/06 22:26:43 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisEvent_h
#define Xyce_N_ANP_AnalysisEvent_h

#include <N_UTL_Dump.h>

namespace Xyce {
namespace Analysis {

struct AnalysisEvent
{
  friend class Dump<AnalysisEvent>;

public:
  enum State {INITIALIZE, STEP_STARTED, STEP_SUCCESSFUL, STEP_FAILED, FINISH};
  enum OutputType {DC, TRAN, AC, AC_IC, HB_FD, HB_TD, HB_IC, HB_STARTUP, DCOP, HOMOTOPY, MPDE, MPDE_IC, SENS};

  AnalysisEvent(State state, OutputType output_type, double step = 0.0, int count = 0)
    : state_(state),
      outputType_(output_type),
      step_(step),
      count_(count)
  {}

  const State           state_;
  const OutputType      outputType_;
  const double          step_;
  const int             count_;
};

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::State &state);

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::OutputType &type);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_AnalysisEvent_h
