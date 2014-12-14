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
// Filename       : $RCSfile: N_IO_OutputResults.C,v $
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
// Revision Number: $Revision: 1.1.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_SweepParam.h>
#include <N_ERH_Message.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputResults.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace IO {

typedef std::list<Util::Param> ParameterList;

//-----------------------------------------------------------------------------
// Class         : OutputMOR
// Purpose       : Output class for OutputMOR runs
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : OutputMOR::OutputMOR
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OutputResults::OutputResults(const std::string &netlist_filename)
  : netlistFilename_(netlist_filename),
    os_(0),
    noIndexResult_(false)
{}

OutputResults::~OutputResults()
{
  delete os_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO OutputMgr::addResultParams
// Purpose       : Sets the RESULT calculation parameters.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/29/04
//-----------------------------------------------------------------------------
bool
OutputResults::addResultParams(
  const Util::OptionBlock &     option_block)
{
// first check to see that there is only 1 PARAM set.  They need to be in
// separate lines for this to work.
  int countPar = 0;
  for (ParameterList::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end(); it != end; ++it)
  {
    if ((*it).uTag() == "EXPRESSION")
      ++countPar;
  }

  if (countPar > 1)
  {
    Report::UserFatal0() << "Only one expression per .RESULT command.  Each parameter needs its own .RESULT line.";
  }

  Util::ExpressionData * expDataPtr;

  for (ParameterList::const_iterator it = option_block.getParams().begin(), end = option_block.getParams().end(); it != end; ++it)
  {
    const Util::Param &expParam = *it;

    if (!expParam.hasExpressionValue())
    {
      Report::DevelFatal0() << "Parameter must be an expression in .RESULT command";
    }
    else
    {
// expression should have already been resolved.  Check for this
// case before we try and create a new expression.
      if (expParam.getType() == Util::EXPR)
      {
        expDataPtr = new Util::ExpressionData(expParam.getValue<Util::Expression>());
      }
      else
      {
        expDataPtr = new Util::ExpressionData(expParam.stringValue());
      }
      resultVector_.push_back(expDataPtr);
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputRESULT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputResults::output(
  Parallel::Machine             comm,
  OutputMgr &                   output_manager,
  double                        circuit_time,
  const Analysis::SweepVector & step_sweep_vector,
  int                           step_loop_number,
  const N_LAS_Vector &          solution_vector,
  const N_LAS_Vector &          state_vector,
  const N_LAS_Vector &          store_vector)
{
  std::string delim = " ";
  int width = 20;
  int precision = 8;

  if (Parallel::rank(comm) == 0)
  {
    if (!os_)
    {
      std::string resultfilename = netlistFilename_ + ".res";

      os_ = new std::ofstream(resultfilename.c_str());

      os_->setf(std::ios::scientific);
      os_->precision(precision);

      if (!noIndexResult_)
        (*os_) << "STEP";

      for (Analysis::SweepVector::const_iterator it = step_sweep_vector.begin(), end = step_sweep_vector.end(); it != end; ++it)
      {
        (*os_) << delim << std::setw(width) << (*it).name;
      }

      for (ResultVector::const_iterator it = resultVector_.begin(), end = resultVector_.end(); it != end; ++it)
      {
        const Util::ExpressionData *expDataPtr = (*it);
        (*os_) << delim << std::setw(width) << expDataPtr->getExpression();
      }
      (*os_) << std::endl;
    }
  }
  
  if (Parallel::rank(comm) == 0)
  {
    os_->setf(std::ios::left, std::ios::adjustfield);
    if (!noIndexResult_)
      (*os_) << step_loop_number;

    for (Analysis::SweepVector::const_iterator it = step_sweep_vector.begin(), end = step_sweep_vector.end(); it != end; ++it)
    {
      (*os_) << delim << std::setw(width) << (*it).currentVal;
    }
  }

  for (ResultVector::const_iterator it = resultVector_.begin(), end = resultVector_.end(); it != end; ++it)
  {
    Util::ExpressionData *expDataPtr = (*it);

    expDataPtr->setup(comm, output_manager);

    double result = expDataPtr->evaluate(comm, circuit_time, &solution_vector, &state_vector, &store_vector);
    if (Parallel::rank(comm) == 0)
    {
      (*os_) << delim << std::setw(width) << result;
    }
  }

  if (Parallel::rank(comm) == 0)
    (*os_) << std::endl;
}


void
OutputResults::steppingComplete()
{
  // Deal with the result file:
  if (os_)
  {
    (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  delete os_;
  os_ = 0;
}

} // namespace IO
} // namespace Xyce
