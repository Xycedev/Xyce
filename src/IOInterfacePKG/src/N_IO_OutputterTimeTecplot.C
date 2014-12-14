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
// Filename       : $RCSfile: N_IO_OutputterTimeTecplot.C,v $
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
// Revision Date  : $Date: 2014/08/26 22:46:45 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterTimeTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : TimeTecplot::TimeTecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeTecplot::TimeTecplot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeTecplot::~TimeTecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeTecplot::~TimeTecplot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeTecplot::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecplot::doOutputTime(
  Parallel::Machine           comm,
  const N_LAS_Vector &        solnVecPtr,
  const N_LAS_Vector &        stateVecPtr,
  const N_LAS_Vector &        storeVecPtr)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    os_->setf(std::ios::scientific);
    os_->precision(printParameters_.streamPrecision_);
    os_->setf(std::ios::left, std::ios::adjustfield);
  }

  if (os_ && index_ == 0)
    tecplotTimeHeader(*os_, currentStep_ == 0, outputManager_.getNetListFilename() + " - " + outputManager_.getTitle(), opList_, outputManager_);

  for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    double result = getValue(comm, *(*it), Util::Op::OpData(0, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0)).real();
    result = filter(result, printParameters_.filter_);
    if ((*it)->id() == Util::Op::identifier<OutputMgrTimeOp>())
      result *= printParameters_.outputTimeScaleFactor_;

    if (os_)
      (*os_) << std::setw(printParameters_.streamWidth_) << result << " ";
  }

  if (os_)
    (*os_) << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : TimeTecplot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecplot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintEndOfSimulationLine())
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeTecplot:::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeTecplot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : TimeTecplot:: doSteppingComplete()
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecplot::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintEndOfSimulationLine())
      (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
