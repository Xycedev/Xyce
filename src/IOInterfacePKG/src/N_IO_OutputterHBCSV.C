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
// Filename       : $RCSfile: N_IO_OutputterHBCSV.C,v $
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
// Revision Number: $Revision: 1.6.2.3 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHBCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {
namespace Outputter {


//-----------------------------------------------------------------------------
// Class         : HBCSV
// Purpose       : Outputter class for HB runs, CSV (comma separated) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HBCSV::HBCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HBCSV::HBCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    freqPrintParameters_(freq_print_parameters),
    timePrintParameters_(time_print_parameters),
    timeFilename_(),
    freqFilename_(),
    tos_(0),
    fos_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.csv";

  if (freqPrintParameters_.defaultExtension_.empty())
    freqPrintParameters_.defaultExtension_ = ".HB.FD.csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), timePrintParameters_, timeOpList_);
  fixupColumns(comm, outputManager_.getOpBuilderManager(), freqPrintParameters_, freqOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::~HBCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HBCSV::~HBCSV()
{
  outputManager_.closeFile(tos_);
  outputManager_.closeFile(fos_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
  deleteList(freqOpList_.begin(), freqOpList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HBCSV::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doOutputHB(
  Parallel::Machine             comm,
  const std::vector<double> &   timePoints,
  const std::vector<double> &   freqPoints,
  const N_LAS_BlockVector &     timeDomainSolutionVec,
  const N_LAS_BlockVector &     freqDomainSolutionVecReal,
  const N_LAS_BlockVector &     freqDomainSolutionVecImaginary,
  const N_LAS_BlockVector &     timeDomainStoreVec,
  const N_LAS_BlockVector &     freqDomainStoreVecReal,
  const N_LAS_BlockVector &     freqDomainStoreVecImaginary)
{
  int index = 0;

  int blockCount = timeDomainSolutionVec.blockCount();

  if (Parallel::rank(comm) == 0 && !fos_ && !tos_)
  {
    timeFilename_ = outputFilename(timePrintParameters_.filename_, timePrintParameters_.defaultExtension_, timePrintParameters_.suffix_, outputManager_.getNetListFilename());
    freqFilename_ = outputFilename(freqPrintParameters_.filename_, freqPrintParameters_.defaultExtension_, freqPrintParameters_.suffix_, outputManager_.getNetListFilename());

    tos_ = outputManager_.openFile(timeFilename_);
    fos_ = outputManager_.openFile(freqFilename_);

    printHeader(*tos_, timePrintParameters_);
    printHeader(*fos_, freqPrintParameters_);
  }

  // Loop over the time points of the N_LAS_BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);
    outputManager_.setCircuitFrequency(freqPoints[iblock]);

    N_LAS_Vector * solutionVecPtr = &(timeDomainSolutionVec.block(iblock));
    N_LAS_Vector * realVecPtr = &(freqDomainSolutionVecReal.block(iblock));
    N_LAS_Vector * imagVecPtr = &(freqDomainSolutionVecImaginary.block(iblock));

    { // periodic time-domain steady-state output
      int column_index = 0;
      for (Util::Op::OpList::const_iterator it = timeOpList_.begin(); it != timeOpList_.end(); ++it, ++column_index)
      {
        double result = getValue(comm, *(*it), Util::Op::OpData(index_, solutionVecPtr, 0, 0, 0, 0)).real();
        if (tos_)
          printValue(*tos_, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);
      }
    }

    { // Fourier coefficient output
      int column_index = 0;
      for (Util::Op::OpList::const_iterator it = freqOpList_.begin(); it != freqOpList_.end(); ++it, ++column_index)
      {
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        double result = getValue(comm, *(*it), Util::Op::OpData(index_, realVecPtr, imagVecPtr, 0, 0, 0)).real();
        if (fos_)
          printValue(*fos_, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, result);
      }
    }

    if (fos_)
      *fos_ << std::endl;
    if (tos_)
      *tos_ << std::endl;

    ++index_;
  }
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doFinishOutput()
{
  if (fos_)
  {
    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(fos_);
      fos_ = 0;
    }
  }
  if (tos_)
  {
    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(tos_);
      tos_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HBCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doSteppingComplete()
{
  outputManager_.closeFile(tos_);
  tos_ = 0;
  outputManager_.closeFile(fos_);
  fos_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
