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
// Filename       : $RCSfile: N_IO_OutputterTimeCSV.C,v $
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
// Revision Number: $Revision: 1.7.2.3 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : TimeCSV
// Purpose       : Outputter class for transient runs and "CSV"
//                 (comma separarated variable) output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeCSV::TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::TimeCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::~TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::~TimeCSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : TimeCSV::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::timeHeader(Parallel::Machine comm)
{
  if (Parallel::rank(comm) == 0)
  {
    printHeader(*os_, printParameters_);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPRINT_
// Purpose       : .PRINT output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void TimeCSV::doOutputTime(
  Parallel::Machine     comm,
  const N_LAS_Vector &  solnVecPtr,
  const N_LAS_Vector &  stateVecPtr,
  const N_LAS_Vector &  storeVecPtr)
{
  double time = outputManager_.getCircuitTime();

  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    timeHeader(comm);
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

void TimeCSV::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::doSteppingComplete()
{
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
