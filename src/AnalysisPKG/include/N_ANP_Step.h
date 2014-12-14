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
// Filename       : $RCSfile: N_ANP_Step.h,v $
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
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/07/16 17:32:40 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Step_h
#define Xyce_N_ANP_Step_h

#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : Step
// Purpose       : Step analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Step : public AnalysisBase
{
public:
  Step( AnalysisManager &analysis_manager, AnalysisBase &main_analysis )
    : AnalysisBase(analysis_manager),
      mainAnalysis_(main_analysis),
      stepLoopSize_(0),
      stepLoopIter_(0),
      stepParamVec_()
  {}

  virtual ~Step()
  {}

  bool setAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

  void setParamsWithOutputMgrAdapter(OutputMgrAdapter & outputManagerAdapter)
  {
    outputManagerAdapter_.setStepParamVec( stepParamVec_ );
  }

  int getStepIter()
  { return stepLoopIter_; }

  virtual bool getDCOPFlag ();

  virtual bool run();
  virtual bool init();
  virtual bool loopProcess();
  virtual bool processSuccessfulStep();
  virtual bool processFailedStep();
  virtual bool finish();
  virtual bool handlePredictor() { return true; }

private:
  AnalysisBase &              mainAnalysis_;
  std::vector <SweepParam>    stepParamVec_;
  int stepLoopSize_;
  int stepLoopIter_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::Step N_ANP_Step;

#endif // Xyce_N_ANP_Step_h
