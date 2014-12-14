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
// Filename      : $RCSfile: N_ANP_Report.C,v $
//
// Purpose       : Contains the class definition for the N_ERH_ErrorMgr
//                 class.  This class is the error handler for Xyce.
//
// Special Notes : 
//
// Creator       : Eric Keiter, SNL,  Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2014/07/16 20:59:26 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <N_ANP_Report.h>
#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {

class messageHeader {
public:
  messageHeader(const AnalysisBase &analysis_base)
    : analysisBase_(analysis_base)
  {}

  const AnalysisBase &          analysisBase_;
};


UserWarning::UserWarning(const AnalysisBase &analysis_base)
  : Report::UserWarning()
{}

UserWarning0::UserWarning0(const AnalysisBase &analysis_base)
  : Report::UserWarning0()
{}

UserFatal::UserFatal(const AnalysisBase &analysis_base)
  : Report::UserFatal()
{}

UserFatal0::UserFatal0(const AnalysisBase &analysis_base)
  : Report::UserFatal0()
{}

DevelFatal::DevelFatal(const AnalysisBase &analysis_base)
  : Report::DevelFatal()
{}

DevelFatal0::DevelFatal0(const AnalysisBase &analysis_base)
  : Report::DevelFatal0()
{}

} // namespace Analysis
} // namespace Xyce
