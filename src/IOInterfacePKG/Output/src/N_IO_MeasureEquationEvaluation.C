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
// Filename      : $RCSfile: N_IO_MeasureEquationEvaluation.C,v $
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

#include <N_IO_MeasureEquationEvaluation.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::EquationEvaluation()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
EquationEvaluation::EquationEvaluation( const Util::OptionBlock & measureBlock):
  Base(measureBlock)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
}

void EquationEvaluation::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();
  
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for relative error measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  outVarValues_.resize( numOutVars_, 0.0 );
 
}

//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void EquationEvaluation::reset() 
{
  initialized_=false;
  calculationDone_=false;
  calculationResult_=0.0;
}


//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void EquationEvaluation:: updateTran(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
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

    // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
    // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
    // take the first element.
    calculationResult_=outVarValues_[0];
  }

}


//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void EquationEvaluation::updateDC(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  // update our outVarValues_ vector
  for( int i=0; i< numOutVars_; i++ )
  {
    outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, stateVec, storeVec, 0);
  }
  // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
  // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
  // take the first element.
  calculationResult_=outVarValues_[0];

}


//-----------------------------------------------------------------------------
// Function      : EquationEvaluation::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2014
//-----------------------------------------------------------------------------
void EquationEvaluation::updateAC(Parallel::Machine comm, const double frequency, const N_LAS_Vector * solnVec, const N_LAS_Vector *imaginaryVec)
{
  if( !calculationDone_ && withinFromToWindow( frequency ) )
  {
    // we're in the frequency window, now we need to calculate the value of this
    // measure and see if it triggers any specified rise, fall, cross windowing.

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0, imaginaryVec);
    }
    // not intuitive, but the output of this measure is just the outputVars_ operator evaluated 
    // within the FromToWindow.  At this time there shouldn't be more than one outVarValues_ so just
    // take the first element.
    calculationResult_=outVarValues_[0];
 
  }


}

} // namespace Measure
} // namespace IO
} // namespace Xyce
