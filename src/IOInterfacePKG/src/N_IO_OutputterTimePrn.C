//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_OutputterTimePrn.C,v $
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
// Revision Number: $Revision: 1.10.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

static Interface *factory(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters) {
  return new TimePrn(comm, output_manager, print_parameters);
}

static void registerOutputter(OutputMgr &output_manager) {
  output_manager.registerOutputter(Analysis::ANP_MODE_DC_SWEEP, OutputType::DC, Format::STD, factory);
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::TimePrn
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
TimePrn::TimePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".prn";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::~TimePrn
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
TimePrn::~TimePrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doOutputTime
// Purpose       : Output the current data at a time point
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void
TimePrn::doOutputTime(
  Parallel::Machine     comm,
  const N_LAS_Vector &  solnVecPtr,
  const N_LAS_Vector &  stateVecPtr,
  const N_LAS_Vector &  storeVecPtr)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    printHeader(*os_, printParameters_);
  }

  int column_index = 0;
  for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(comm, *(*it), Util::Op::OpData(index_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0)).real();
    result = filter(result, printParameters_.filter_);
    if ((*it)->id() == Util::Op::identifier<OutputMgrTimeOp>())
      result *= printParameters_.outputTimeScaleFactor_;

    if (os_)
      printValue(*os_, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);
  }

  if (os_)
    *os_ << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doFinishOutput
// Purpose       : Output the footer, close stream
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintEndOfSimulationLine())
        *os_ << "End of Xyce(TM) Simulation" << std::endl;

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimePrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doSteppingComplete
// Purpose       : output footer and close stream for parameter sweep file
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintEndOfSimulationLine())
      *os_ << "End of Xyce(TM) Parameter Sweep" << std::endl;

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
