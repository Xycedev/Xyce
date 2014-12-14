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
// Filename       : $RCSfile: N_IO_OutputterSensitivityTecplot.C,v $
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
// Revision Number: $Revision: 1.9.2.4 $
//
// Revision Date  : $Date: 2014/09/02 22:39:44 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterSensitivityTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : SensitivityTecPlot
// Purpose       : Outputter class for sensitivity output, TecPlot output
//                 format
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
SensitivityTecPlot::SensitivityTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = "SENS.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::~SensitivityTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
SensitivityTecPlot::~SensitivityTecPlot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doOutputSensitivity
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void
SensitivityTecPlot::doOutputSensitivity(
  Parallel::Machine             comm,
  const std::vector<double> &   objective_values,
  const std::vector<double> &   direct_values,
  const std::vector<double> &   adjoint_values,
  const std::vector<double> &   scaled_direct_values,
  const std::vector<double> &   scaled_adjoint_values,
  const N_LAS_Vector &          solution_vector,
  const N_LAS_Vector &          state_vector,
  const N_LAS_Vector &          store_vector)
{
  double tmpTime = outputManager_.getCircuitTime(); // outputManager_.getAnaIntPtr()->getTime();

  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);

    (*os_).setf(std::ios::left, std::ios::adjustfield);
  }

  if (os_ && index_ == 0)
    tecplotTimeHeader(*os_, currentStep_ == 0, outputManager_.getNetListFilename() + " - " + outputManager_.getTitle(), opList_, outputManager_);

  int column_index = 0;
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(comm, *(*it),
             Util::Op::OpData(index_, &solution_vector, 0, &state_vector, &store_vector, 0, 
             &objective_values, &direct_values, &scaled_direct_values, &adjoint_values, &scaled_adjoint_values))
                    .real();

    result = filter(result, printParameters_.filter_);

    if ((*it)->id() == Util::Op::identifier<OutputMgrTimeOp>())
      result *= printParameters_.outputTimeScaleFactor_;

    if (Parallel::rank(comm) == 0)
      printValue(*os_, printParameters_.table_.columnList_[column_index],
                 printParameters_.delimiter_, column_index, result);
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doFinishOutput()
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
// Function      : SensitivityTecPlot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
SensitivityTecPlot::doStartStep(
  int                           step,
  int                           max_step)
{
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doSteppingComplete()
{
  // close the sensitivity file.
  if (os_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*os_) << "End of Xyce(TM) Sensitivity Simulation" << std::endl;
    }
  }

  outputManager_.closeFile(os_);
  os_ = 0;
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce
