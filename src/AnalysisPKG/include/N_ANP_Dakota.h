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
// Filename       : $RCSfile: N_ANP_Dakota.h,v $
//
// Purpose        : Dakota analysis class
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
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/07/16 17:32:40 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Dakota_h
#define Xyce_N_ANP_Dakota_h

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {
//-------------------------------------------------------------------------
// Class         : Dakota
// Purpose       : Dakota analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Dakota : public AnalysisBase
{
  public:
    Dakota(AnalysisManager &anaManagerPtr, AnalysisBase &anaType );

    virtual ~Dakota( );
  

    virtual bool getDCOPFlag () {return true;}// this needs to be fixed.
    virtual bool run();
    virtual bool init() { return true; }
    virtual bool loopProcess() { return true; }
    virtual bool processSuccessfulStep() { return true; }
    virtual bool processFailedStep() { return true; }
    virtual bool finish() { return true; }
    virtual bool handlePredictor() { return true; }

private:
    AnalysisBase &              mainAnalysis_;

  // these should be part of this class, but are still too entangled in
  // the AnalysisManager class to be moved
  //vector <N_TIA_SweepParam> stepParamVec_;
  //
  //bool stepLoopInitialized_;    // true if the step loop has been set up.
  //int stepLoopSize_;
  //int stepLoopIter_;

};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::Dakota N_ANP_Dakota;

#endif // Xyce_N_ANP_Dakota_h
