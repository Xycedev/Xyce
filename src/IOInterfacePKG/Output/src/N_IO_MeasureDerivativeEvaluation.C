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
// Filename      : $RCSfile: N_IO_MeasureDerivativeEvaluation.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.14.2.1 $
// Revision Date  : $Date: 2014/08/28 22:41:46 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureDerivativeEvaluation.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::DerivativeEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
DerivativeEvaluation::DerivativeEvaluation( const Util::OptionBlock & measureBlock):
  Base(measureBlock),
  firstTimeValue_(0.0),
  lastTimeValue_(0.0),
  firstSignalValue_(0.0),
  lastSignalValue_(0.0)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
}

void DerivativeEvaluation::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for statistical measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void DerivativeEvaluation::reset() 
{
  initialized_=false;
  calculationDone_=false;
  firstTimeValue_ = 0.0;
  firstSignalValue_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateTran(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  if( !calculationDone_ && withinFromToWindow( circuitTime ) )
  {
    // we're in the transient window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, stateVec, storeVec, 0);
    }

    if( !initialized_  )
    {
      // Calculating an average slope by calculating and averaging the slope
      // between successive paris of points
      firstTimeValue_ = circuitTime;
      firstSignalValue_ = outVarValues_[0];
      initialized_ = true;
    }

    lastTimeValue_ = circuitTime;
    lastSignalValue_ = outVarValues_[0];
  }
}


//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void DerivativeEvaluation::updateDC(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{

}

//-----------------------------------------------------------------------------
// Function      : DerivativeEvaluation::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double DerivativeEvaluation::getMeasureResult()
{
  if( initialized_ )
  {
    calculationResult_ =  (lastSignalValue_ - firstSignalValue_) / (lastTimeValue_ - firstTimeValue_);
  }
  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
