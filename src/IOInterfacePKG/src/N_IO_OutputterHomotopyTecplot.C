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
// Filename       : $RCSfile: N_IO_OutputterHomotopyTecplot.C,v $
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
// Revision Number: $Revision: 1.8.2.3 $
//
// Revision Date  : $Date: 2014/09/02 22:39:44 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHomotopyTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : HomotopyTecPlot
// Purpose       : Outputter class for homotopy output, TecPlot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyTecPlot::HomotopyTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
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
// Function      : HomotopyTecPlot::~HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyTecPlot::~HomotopyTecPlot()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyTecPlot::homotopyHeader(
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

  std::ostream &os = *os_;

  index_ = 0;
  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();

  if (currentStep_ == 0)
  {
    os << " TITLE = \" Xyce homotopy data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the continuation parameters:
    std::vector<std::string>::const_iterator iter_name;
    for (iter_name = parameter_names.begin(); iter_name!= parameter_names.end(); ++iter_name)
    {
      os << "\" ";
      os << *iter_name;
      os << "\" " << std::endl;
    }

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
    {
      os << "\" " << (*it)->getName() << "\" " << std::endl;
    }
  }

  // output some AUXDATA
  os << "DATASETAUXDATA ";
  os << getTecplotTimeDateStamp();
  os << std::endl;

  os << "ZONE F=POINT";


  if (outputManager_.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin();
         it != outputManager_.getStepParamVec().end(); ++it)
    {
      static const int tecplotHeaderPrecision = 2;
      os.setf(std::ios::scientific);
      os.precision(tecplotHeaderPrecision);
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyTecPlot::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const N_LAS_Vector &                  solution_vector)
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
    homotopyHeader(parameter_names, parameter_values, solution_vector);

  int column_index = 0;
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(comm, *(*it), Util::Op::OpData(0, &solution_vector, 0, 0, 0, 0)).real();
    if (os_)
    {
      if (column_index == 0)
      {
        for (int i = 0; i < parameter_values.size(); ++i)
        {
          printValue(*os_, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);
        }
      }

      printValue(*os_, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doFinishOutput()
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
// Function      : HomotopyTecPlot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyTecPlot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintEndOfSimulationLine() )
      (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
