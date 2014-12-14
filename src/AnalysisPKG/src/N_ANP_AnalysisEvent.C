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
// Filename      : $RCSfile: N_ANP_AnalysisEvent.C,v $
// Purpose       : AnalysisEvent info
// Special Notes :
// Creator       :
// Creation Date :
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.1 $
// Revision Date  : $Date: 2014/08/06 22:26:43 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <ostream>

#include <N_ANP_AnalysisEvent.h>

namespace Xyce {
namespace Analysis {

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::State &state)
{
  switch (state) {
    case AnalysisEvent::INITIALIZE:
      os << "INITIALIZE";
      break;

    case AnalysisEvent::STEP_STARTED:
      os << "STEP_STARTED";
      break;

    case AnalysisEvent::STEP_SUCCESSFUL:
      os << "STEP_SUCCESSFUL";
      break;

    case AnalysisEvent::STEP_FAILED:
      os << "STEP_FAILED";
      break;

    case AnalysisEvent::FINISH:
      os << "FINISH";
      break;
  }

  return os;
}

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::OutputType &type)
{
  switch (type) {

    case AnalysisEvent::DC:
      os << "DC";
      break;

    case AnalysisEvent::TRAN:
      os << "TRAN";
      break;

    case AnalysisEvent::AC:
      os << "AC";
      break;

    case AnalysisEvent::AC_IC:
      os << "AC_IC";
      break;

    case AnalysisEvent::HB_FD:
      os << "HB_FD";
      break;

    case AnalysisEvent::HB_TD:
      os << "HB_TD";
      break;

    case AnalysisEvent::HB_IC:
      os << "HB_IC";
      break;

    case AnalysisEvent::HB_STARTUP:
      os << "HB_STARTUP";
      break;

    case AnalysisEvent::DCOP:
      os << "DCOP";
      break;

    case AnalysisEvent::HOMOTOPY:
      os << "HOMOTOPY";
      break;

    case AnalysisEvent::MPDE:
      os << "MPDE";
      break;

    case AnalysisEvent::MPDE_IC:
      os << "MPDE_IC";
      break;

    case AnalysisEvent::SENS:
      os << "SENS";
      break;
  }

  return os;
}

} // namespace Analysis

template<>
std::ostream &Dump<Analysis::AnalysisEvent>::dump(std::ostream &os) const
{
  // os << "AnalysisEvent "
  //    << "state_ " << t_.state_
  //    << ", outputType_ " << t_.outputType_
  //    << ", step_ " << t_.step_
  //    << ", count_ " << t_.count_ << std::endl;

  return os;
}

} // namespace Xyce
