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
// Filename      : $RCSfile: N_IO_MeasureIntegralEvaluation.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.15.2.1 $
// Revision Date  : $Date: 2014/08/28 22:41:46 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureIntegralEvaluation.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::IntegralEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
IntegralEvaluation::IntegralEvaluation( const Util::OptionBlock & measureBlock):
  Base(measureBlock),
  integralValue_(0.0),
  lastTimeValue_(0.0),
  lastSignalValue_(0.0),
  totalAveragingWindow_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
}


void IntegralEvaluation::prepareOutputVariables()
{
  outVarValues_.resize( outputVars_.size(), 0.0 );
  
  // this measurement should have only one dependent variable.
  if (outVarValues_.size() > 1 )
    Xyce::Report::UserError0() << "Too many dependent variables for statistical measure, \"" << name_ << "\"";
}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void IntegralEvaluation::reset() 
{
  initialized_=false;
  calculationDone_=false;
  integralValue_ = 0.0;
  totalAveragingWindow_ =0.0;
}

//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateTran(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  
  if( !calculationDone_ && withinFromToWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    for( int i = 0; i < outVarValues_.size(); i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, stateVec, storeVec, 0);
    }

    if( initialized_  )
    {
      // update integral of outVarValues_[0];
      integralValue_ += 0.5 * (circuitTime - lastTimeValue_) * (outVarValues_[0] + lastSignalValue_);
      totalAveragingWindow_ += (circuitTime - lastTimeValue_);
    }

    lastTimeValue_ = circuitTime;
    lastSignalValue_ = outVarValues_[0];
    initialized_=true;
  }
}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::updateMeasures()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void IntegralEvaluation::updateDC(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{

}


//-----------------------------------------------------------------------------
// Function      : IntegralEvaluation::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double IntegralEvaluation::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  integralValue_;
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
