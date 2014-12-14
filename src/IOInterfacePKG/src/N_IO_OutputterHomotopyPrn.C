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
// Filename       : $RCSfile: N_IO_OutputterHomotopyPrn.C,v $
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
// Revision Number: $Revision: 1.7.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHomotopyPrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {



//-----------------------------------------------------------------------------
// Class         : HomotopyPrn
// Purpose       : Outputter class for homotopy output, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyPrn::HomotopyPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::~HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyPrn::~HomotopyPrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::homotopyHeader(
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const N_LAS_Vector &                  solution_vector)
{
  if (columnList_.empty())
  {
    Table::Justification justification = printParameters_.delimiter_.empty() ?
      Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin();
        it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific,
            printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
    }
  }

  index_ = 0;

  int homotopyParamStartIndex=1;
  if (!printParameters_.printIndexColumn_) // if noindex, then use 0 for start of homotopy params, otherwise 1
  {
    homotopyParamStartIndex=0;
  }

  if (currentStep_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin();
        it != printParameters_.table_.columnList_.end(); ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }

      if (column_index == homotopyParamStartIndex)
      {
        for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
        {
          if (it2 != columnList_.begin())
          {
            *os_ << printParameters_.delimiter_;
          }
          printHeader(*os_, (*it2));
        }
      }

      printHeader(*os_, (*it));
    }

    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const N_LAS_Vector &                  solution_vector)
{
  double tmpTime = outputManager_.getCircuitTime(); // outputManager_.getAnaIntPtr()->getTime();

  if (Parallel::rank(comm) == 0 && !os_)
    {
      outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
      os_ = outputManager_.openFile(outFilename_);

      homotopyHeader(parameter_names, parameter_values, solution_vector);
  }

  int homotopyParamStartIndex=1;
  if (!printParameters_.printIndexColumn_) // if noindex, then use 0 for start of homotopy params, otherwise 1
  {
    homotopyParamStartIndex=0;
  }

  int column_index = 0;
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(comm, *(*it), Util::Op::OpData(index_, &solution_vector, 0, 0, 0, 0)).real();
    if (Parallel::rank(comm) == 0)
    {
      if (column_index == homotopyParamStartIndex)
      {
        for (int i = 0; i < parameter_values.size(); ++i)
        {
          printValue(*os_, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);
        }
      }

      printValue(*os_, printParameters_.table_.columnList_[column_index],
          printParameters_.delimiter_, column_index, result);
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintEndOfSimulationLine())
        (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doSteppingComplete()
{
  // close the homotopy file.
  if (os_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
      (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce
