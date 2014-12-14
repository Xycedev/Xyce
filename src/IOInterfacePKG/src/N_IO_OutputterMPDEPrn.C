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
// Filename       : $RCSfile: N_IO_OutputterMPDEPrn.C,v $
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
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterMPDEPrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>
#include <N_MPDE_Manager.h>

namespace Xyce {
namespace IO {
namespace Outputter {


//-----------------------------------------------------------------------------
// Class         : MPDEPrn
// Purpose       : Outputter class for MPDE runs, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDEPrn::MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDEPrn::MPDEPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    n1_(0),
    n2_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".MPDE.prn";

  printParameters_.table_.addColumn("TIME1", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);
  printParameters_.table_.addColumn("TIME2", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::~MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDEPrn::~MPDEPrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::mpdeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::mpdeHeader()
{}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doOutputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points, 
  const N_LAS_Vector &          solution_vector)
{
  const N_LAS_BlockVector & blockSolVecPtr = dynamic_cast<const N_LAS_BlockVector &>(solution_vector);

  int blockCount = blockSolVecPtr.blockCount();
  n2_ = blockCount; // fast time points.
  ++n1_;            // slow time points.  increments by one each call.

  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    mpdeHeader();
  }

  // Loop over the fast time points of the N_LAS_BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    const N_LAS_Vector &solnVecPtr = (iblock == n2_ ? blockSolVecPtr.block(0) : blockSolVecPtr.block(iblock));

    if (os_)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = time;
      double second = fast_time_points[iblock]; //genie 121913. This should be a bug since fast_time_points has size = n2_ so its index is at most n2_-1

      // time 1:
      printValue(*os_, printParameters_.table_.columnList_[0], printParameters_.delimiter_, 0, first);

      // time 2:
      printValue(*os_, printParameters_.table_.columnList_[1], printParameters_.delimiter_, 1, second);
    }

    int column_index = 2;
    for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
    {
      double result = getValue(comm, *(*it), Util::Op::OpData(0, &solnVecPtr, 0, 0, 0, 0)).real();
      if (os_)
      {
        printValue(*os_, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);
      }
    }

    if (os_)
    {
      (*os_) << std::endl;
    }
  } // fast time scale loop.
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doFinishOutput()
{
  outputManager_.closeFile(os_);
  os_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
MPDEPrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doSteppingComplete()
{}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
