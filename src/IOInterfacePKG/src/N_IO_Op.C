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
// Filename       : $RCSfile: N_IO_Op.C,v $
//
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.31 $
//
// Revision Date  : $Date: 2014/08/11 18:03:53 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_fwd.h>
#include <N_IO_Op.h>
#include <N_DEV_Op.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_IO_MeasureBase.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Marshal.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : CurrentIndexOp::get
// Purpose       : get the "index" from the currently active outputter
// Special Notes : the index is basically the line number of output, starting
//                 at zero for the first line
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
CurrentIndexOp::get(const CurrentIndexOp &op, const Util::Op::OpData &op_data)
{
  return op_data.currentIndex_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrTimeOp::get
// Purpose       : get the current simulation time being output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrTimeOp::get(const OutputMgrTimeOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getTime();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrFrequencyOp::get
// Purpose       : get the current frequency being output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrFrequencyOp::get(const OutputMgrFrequencyOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getFrequency();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrTemperatureOp::get
// Purpose       : get the current simulation temperature
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrTemperatureOp::get(const OutputMgrTemperatureOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getTemperature();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrStepSweepOp::get
// Purpose       : get the current value of the step parameter being swept.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrStepSweepOp::get(const OutputMgrStepSweepOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getStepSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrDCSweepOp::get
// Purpose       : get the current value of the DC voltage being swept.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrDCSweepOp::get(const OutputMgrDCSweepOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getDCSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrDCSweepCurrentValueOp::get
// Purpose       : get the current value of the DC voltage being swept.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
OutputMgrDCSweepCurrentValueOp::get(const OutputMgrDCSweepCurrentValueOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getPRINTDCvalue();
}

//-----------------------------------------------------------------------------
// Function      : ObjectiveOp::get
// Purpose       : get the current value of the objective function
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ObjectiveOp::get(const ObjectiveOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.objective_.var1.empty() && op.objective_.var2.empty())
  {
    // saving single values(potentially all values) as there isn't
    // any external data to tell us what needs to be saved
    result = complex(op.objective_.save(op.comm_, op.outputMgr_.getCircuitTime(), op_data.realSolutionVector_, op_data.stateVector_, op_data.realStoreVector_), 0.0);
  }
  else
  {
    // use user supplied external data to save only the simulation
    // results we really need.
    double v1 = 0.0;
    double v2 = 0.0;
    // Util::Param param;
    // if (!objective.var1.empty())
    // {
    //   param.setTag(objective.var1);
    //   v1 = getPrintValue(outputMgr_, param, realSolutionVector_, stateVector_, storeVector_);
    // }
    // if (!objective.var2.empty())
    // {
    //   param.setTag(objective.var2);
    //   v2 = getPrintValue(outputMgr_, param, realSolutionVector_);
    // }
    result = complex(op.objective_.save(op.comm_, op.outputMgr_.getCircuitTime(), v1, v2, op_data.realSolutionVector_, op_data.stateVector_, op_data.realStoreVector_), 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionOp::get
// Purpose       : get the current value of a solution vector element
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionOp::get(const SolutionOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::get
// Purpose       : get a solution variable in preparation for computing real
//                 part
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionRealOp::get(const SolutionRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::eval
// Purpose       : take the real part of a complex number
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::get
// Purpose       : get a solution variable in preparation for computing
//                 imaginary part
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionImaginaryOp::get(const SolutionImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::get
// Purpose       : get the magnitude of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionMagnitudeOp::get(const SolutionMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::eval
// Purpose       : take the magnitude of a solution vector element
// Special Notes : Actually just takes the magnitude of a given complex
//                 value, does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : SolutionPhaseOp::get
// Purpose       : get the phase of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the phase.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionPhaseOp::get(const SolutionPhaseOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionPhaseOp::eval
// Purpose       : compute the phase of a solution vector element
// Special Notes : Actually just computes the phase of a given complex
//                 value, does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionPhaseOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::get
// Purpose       : get the magnitude (in dB) of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 find the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionDecibelsOp::get(const SolutionDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a solution vector element
// Special Notes : Actually just computes the magnitude of a given complex
//                 value, does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SolutionDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceOp::get(const VoltageDifferenceOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : Computes the difference, but does not take the real part.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceRealOp::get(const VoltageDifferenceRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.realSolutionVector_ == 0 ? 0.0 : (*op_data.realSolutionVector_)[op.index1_]);
  }

  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.realSolutionVector_ == 0 ? 0.0 : (*op_data.realSolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::eval
// Purpose       : take the real part of a voltage difference
// Special Notes : must "get" the difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : Computes the difference, but does not take the imaginary
//                 part.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceImaginaryOp::get(const VoltageDifferenceImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
  }

  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::get
// Purpose       : get the imaginary part of a voltage difference
// Special Notes : Must "get" the difference first.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceMagnitudeOp::get(const VoltageDifferenceMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::eval
// Purpose       : Compute the magnitude of a voltage difference
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferencePhaseOp::get(const VoltageDifferencePhaseOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseOp::eval
// Purpose       : Compute the phase of a voltage difference
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferencePhaseOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceDecibelsOp::get(const VoltageDifferenceDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1)
  {
    result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
  }
  if (op.index2_ != -1)
  {
    result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::eval
// Purpose       : Compute the magnitude of a voltage difference (in dB)
// Special Notes : Must "get" difference first
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
VoltageDifferenceDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(std::abs(result));
}

//-----------------------------------------------------------------------------
// Function      : StateOp::get
// Purpose       : Get a value out of the state vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StateOp::get(const StateOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.stateVector_ == 0 ? 0.0 : (*op_data.stateVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreOp::get
// Purpose       : Get a value out of the Store vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreOp::get(const StoreOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realStoreVector_ == 0 ? 0.0 : (*op_data.realStoreVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreRealOp::get
// Purpose       : get a store variable in preparation for computing real
//                 part
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreRealOp::get(const StoreRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreRealOp::eval
// Purpose       : take the real part of a complex number
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : StoreImaginaryOp::get
// Purpose       : get a store variable in preparation for computing
//                 imaginary part
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreImaginaryOp::get(const StoreImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : StoreMagnitudeOp::get
// Purpose       : get the magnitude of a store vector element
// Special Notes : Actually just gets the store variable.  Does not
//                 take the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreMagnitudeOp::get(const StoreMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreMagnitudeOp::eval
// Purpose       : take the magnitude of a store vector element
// Special Notes : Actually just takes the magnitude of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : StorePhaseOp::get
// Purpose       : get the phase of a store vector element
// Special Notes : Actually just gets the store variable.  Does not
//                 take the phase.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StorePhaseOp::get(const StorePhaseOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StorePhaseOp::eval
// Purpose       : compute the phase of a store vector element
// Special Notes : Actually just computes the phase of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StorePhaseOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : StoreDecibelsOp::get
// Purpose       : get the magnitude (in dB) of a store vector element
// Special Notes : Actually just gets the store variable.  Does not
//                 find the magnitude.  eval does that.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreDecibelsOp::get(const StoreDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a store vector element
// Special Notes : Actually just computes the magnitude of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
StoreDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

//-----------------------------------------------------------------------------
// Function      : SensitivityObjFunctionOp::get
// Purpose       : Get a value out of the SensitivityObjFunction vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SensitivityObjFunctionOp::get(const SensitivityObjFunctionOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.objectiveVector_ == 0 || op.index_ >= op_data.objectiveVector_->size() ? 0.0 : (*op_data.objectiveVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpDirectOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SensitivitydOdpDirectOp::get(const SensitivitydOdpDirectOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpDirectVector_ == 0 || op.index_ >= op_data.dOdpDirectVector_->size() ? 0.0 : (*op_data.dOdpDirectVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpDirectScaledOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SensitivitydOdpDirectScaledOp::get(const SensitivitydOdpDirectScaledOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpDirectScaledVector_ == 0 || op.index_ >= op_data.dOdpDirectScaledVector_->size() ? 0.0 : (*op_data.dOdpDirectScaledVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpAdjointOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SensitivitydOdpAdjointOp::get(const SensitivitydOdpAdjointOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpAdjointVector_ == 0 || op.index_ >= op_data.dOdpAdjointVector_->size() ? 0.0 : (*op_data.dOdpAdjointVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpAdjointScaledOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
SensitivitydOdpAdjointScaledOp::get(const SensitivitydOdpAdjointScaledOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpAdjointScaledVector_ == 0 || op.index_ >= op_data.dOdpAdjointScaledVector_->size() ? 0.0 : (*op_data.dOdpAdjointScaledVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : MeasureOp::get
// Purpose       : Get a .measure result
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
MeasureOp::get(const MeasureOp &op, const Util::Op::OpData &op_data)
{
  complex result(const_cast<Measure::Base &>(op.measure_).getMeasureResult(), 0.0);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::ExpressionOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes pre-constructed expression as second argument
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionOp::ExpressionOp(const std::string &name, Util::Expression &expression, Parallel::Machine comm, const OutputMgr &output_manager)
  : Base(name),
    expressionData_(expression),
    comm_(comm),
    outputMgr_(output_manager)
{
  init(comm, output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::ExpressionOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes string as second argument.  expressionData_
//                 constructor will process the string into an expression
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionOp::ExpressionOp(const std::string &name, const std::string &expression, Parallel::Machine comm, const OutputMgr &output_manager)
  : Base(name),
    expressionData_(expression),
    comm_(comm),
    outputMgr_(output_manager)
{
  init(comm, output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::init
// Purpose       : initialize an expression
// Special Notes : runs the ExpressionData::setup method to resolve
//                 symbols
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void ExpressionOp::init(Parallel::Machine comm, const OutputMgr &output_manager)
{
  expressionData_.setup(comm, output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::get
// Purpose       : evaluate an expression
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex
ExpressionOp::get(const ExpressionOp &op, const Util::Op::OpData &op_data)
{
  complex result(op.expressionData_.evaluate(op.comm_, op.outputMgr_.getCircuitTime(), op_data.realSolutionVector_, op_data.stateVector_, op_data.realStoreVector_, op_data.imaginarySolutionVector_), 0.0);

  return result;
}

namespace {

NodeNamePairMap::const_iterator findNode(const std::string &name, const NodeNamePairMap &node_map, const AliasNodeMap &alias_map)
{
  NodeNamePairMap::const_iterator node_it = node_map.find(name);
  if (node_it == node_map.end()) {
    AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
    if (alias_node_it != alias_map.end())
      node_it = node_map.find((*alias_node_it).second);
  }

  return node_it;
}

} // namespace <unnamed>

} // namespace IO
} // namespace Xyce
