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
// Filename       : $RCSfile: N_IO_OutputterTimeProbe.C,v $
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
// Revision Number: $Revision: 1.9.2.3 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>
#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace IO {
namespace Outputter {


//-----------------------------------------------------------------------------
// Class         : TimeProbe
// Purpose       : Outputter class for transient runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeProbe::TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeProbe::TimeProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    printCount_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csd";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::~TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeProbe::~TimeProbe()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::timeHeader(Parallel::Machine comm)
{
  if (os_)
  {
    std::ostream &os = *os_;

    // count the number of output variables.
    printCount_ = 0;
    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
    {
      if ((*it)->id() != Util::Op::identifier<Util::Op::UndefinedOp>())
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << Util::Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;

    os.setf(std::ios::scientific);
    os.precision(0); // streamPrecision_);
    if (outputManager_.getStepParamVec().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      os << "SUBTITLE='Step param";
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
    }

    os << " ' " << std::endl;;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
      os << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
    else
      os << "ANALYSIS='DC Sweep' " << "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='NO' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      os << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" <<
        std::endl;
    }
    else
    {
      os << "SWEEPVAR='";
      os << outputManager_.getPRINTDCname();
      os << "' SWEEPMODE=";
      if ((outputManager_.getDCParamVec().size() > 0) && (outputManager_.getDCParamVec()[0].type == "LIST"))
      {
        os << "'LIST'" << std::endl;
      }
      else
      {
        os << "'VAR_STEP'" << std::endl;
      }
    }

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      os << "XBEGIN='" << outputManager_.getAnaIntPtr()->getInitialTime()
         << "' XEND='" << outputManager_.getAnaIntPtr()->getFinalTime() << "'" << std::endl;
    }
    else
    {
      os << "XBEGIN='" << outputManager_.getPRINTDCstart()
         << "' XEND='" << outputManager_.getPRINTDCstop() << "'" << std::endl;
    }

    os << "FORMAT='0 VOLTSorAMPS;EFLOAT : "
       << "NODEorBRANCH;NODE  '  " << std::endl;
    os << "DGTLDATA='NO'";

    int dcSize = outputManager_.getDCParamVec().size();
    if (dcSize > 1)
    {
      os << "  ";
      for (int idc=1;idc<dcSize;++idc)
      {
        os << "SWEEP" << idc+1 << "PARM='";
        os << outputManager_.getDCParamVec()[idc].name;
        os << "' ";
        os << "SWEEP" << idc+1 << "VALUE='";
        os.setf(std::ios::scientific);
        os.precision(printParameters_.streamPrecision_);
        os << outputManager_.getDCParamVec()[idc].currentVal;
        os << "' ";
        os << std::endl;
      }
    }
    else
    {
      os << std::endl;
    }
    os << "#N" << std::endl;

    int i = 0;
    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      if (i > 18)
      {
        i = 0;
        os << std::endl;
      }
      os << "'" << (*it)->getName() << "' ";
    }
    if (i != 0)
      os << std::endl;

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(std::ios::left, std::ios::adjustfield);

  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeProbe::doOutputTime(
  Parallel::Machine     comm,
  const N_LAS_Vector &  solnVecPtr,
  const N_LAS_Vector &  stateVecPtr,
  const N_LAS_Vector &  storeVecPtr)
{
  if (Parallel::rank(comm) == 0 && !os_)
    {
      outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
      os_ = outputManager_.openFile(outFilename_);
      timeHeader(comm);
    }

  if (os_) {
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
      *os_ << "#C " << outputManager_.getCircuitTime() << " " << printCount_ << std::endl;
    else
      *os_ << "#C " << outputManager_.getPRINTDCvalue() << " " << printCount_ << std::endl;
  }

  int i = 1;
  for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
  {
    double result = getValue(comm, *(*it), Util::Op::OpData(index_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0)).real();
    if ((*it)->id() == Util::Op::identifier<OutputMgrTimeOp>())
      result *= printParameters_.outputTimeScaleFactor_;

    if (os_)
      *os_ << result << ":" << i << (i%5 == 0 ? "\n" : "   ");
  }

  if (os_ && i%5 != 0)
    *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doFinishOutput()
{
  if (os_)
  {
    (*os_) << "#;" << std::endl;

    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeProbe::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doSteppingComplete()
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
