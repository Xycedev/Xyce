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
// Filename       : $RCSfile: N_IO_OutputterMPDETecplot.C,v $
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
// Revision Date  : $Date: 2014/09/02 22:39:44 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterMPDETecplot.h>
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
// Class         : MPDETecplot
// Purpose       : Outputter class for MPDE runs, Tecplot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDETecplot::MPDETecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
MPDETecplot::MPDETecplot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    n1_(0),
    n2_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".MPDE.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::~MPDETecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDETecplot::~MPDETecplot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::mpdeHeader()
{
  (*os_) << " TITLE = \" Xyce MPDE data, " << outputManager_.getNetListFilename() << "\", " << std::endl
                   << "\tVARIABLES = \"T1(sec) \", \"T2(sec)\", " << std::endl;

  // output the user-specified solution vars:
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
  {
    (*os_) << "\" "<< (*it)->getName() << "\" " << std::endl;
  }

  // output some AUXDATA
  (*os_) << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl
                   << "ZONE I=" << n2_ + 1 << ", " << " J=" << n1_ << ", " << " F=POINT\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
MPDETecplot::doOutputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points, 
  const N_LAS_Vector &          solution_vector)
{
  const N_LAS_BlockVector &blockSolVecPtr = dynamic_cast<const N_LAS_BlockVector &>(solution_vector);

  int blockCount = blockSolVecPtr.blockCount();
  n2_ = blockCount; // fast time points.
  ++n1_;            // slow time points.  increments by one each call.

  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    os_->setf(std::ios::scientific);
    os_->precision(printParameters_.streamPrecision_);
    os_->setf(std::ios::left, std::ios::adjustfield);
  }

  if (os_ && index_ == 0)
    mpdeHeader();

  // Loop over the fast time points of the N_LAS_BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    const N_LAS_Vector &solnVecPtr = (iblock == n2_ ? blockSolVecPtr.block(0) : blockSolVecPtr.block(iblock));

    if (Parallel::rank(comm) == 0)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = 0.0;
      double second = 0.0;

      second = fast_time_points[iblock];
      first  = time;

      // time 1:
      (*os_) << std::setw(printParameters_.streamWidth_) << first << " " << std::setw(printParameters_.streamWidth_) << second;
    }

    for (Util::Op::OpList::const_iterator it = opList_.begin(), end = opList_.end(); it != end; ++it)
    {
      double result = getValue(comm, *(*it), Util::Op::OpData(0, &solnVecPtr, 0, 0, 0, 0)).real();
      if (os_)
        (*os_) << result;
    }

    if (os_)
      (*os_) << "\n";
  }

  if (os_)
    (*os_) << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doFinishOutput()
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
// Function      : MPDETecplot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
MPDETecplot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doSteppingComplete()
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
