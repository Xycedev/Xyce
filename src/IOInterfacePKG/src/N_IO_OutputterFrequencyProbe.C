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
// Filename       : $RCSfile: N_IO_OutputterFrequencyProbe.C,v $
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

#include <N_IO_OutputterFrequencyProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {
//-----------------------------------------------------------------------------
// Class         : FrequencyProbe
// Purpose       : Outputter class for frequency-domain runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyProbe::FrequencyProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
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
// Function      : FrequencyProbe::~FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyProbe::~FrequencyProbe()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::frequencyHeader(Parallel::Machine comm)
{
  if (os_) {
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

    os << "ANALYSIS='AC Sweep' " <<
      "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='YES' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    os << "SWEEPVAR='";
    std::string varName = outputManager_.getPRINTDCname();
    if( varName == "" )
    {
      varName="FREQ";
    }
    os << varName;
    os << "' SWEEPMODE=";
    if ( (!outputManager_.getDCParamVec().empty()) && (outputManager_.getDCParamVec()[0].type == "LIST") )
    {
      os << "'LIST'" << std::endl;
    }
    else
    {
      os << "'VAR_STEP'" << std::endl;
    }


    os << "XBEGIN='" << outputManager_.getPRINTDCstart()
       << "' XEND='" << outputManager_.getPRINTDCstop() << "'" << std::endl;

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
      os << "'" << (*it)->getName() << "' ";
      if (i > 3)
      {
        i = 0;
        os << std::endl;
      }
    }

    if (i != 0)
      os << std::endl;

    os.flush();

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(std::ios::left, std::ios::adjustfield);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyProbe::doOutputFrequency(
  Parallel::Machine     comm,
  double                frequency,
  const N_LAS_Vector &  real_solution_vector,
  const N_LAS_Vector &  imaginary_solution_vector)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openFile(outFilename_);
    frequencyHeader(comm);
  }

  if (os_)
    *os_ << "#C " << frequency << " " << printCount_ << std::endl;

  int i = 1;
  for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
  {
    complex result = getValue(comm, *(*it), Util::Op::OpData(0, &real_solution_vector, &imaginary_solution_vector, 0, 0, 0));
    if (os_)
    {
      *os_ << result.real() << "/" << result.imag() << ":" << i << (i%5 == 0 ? "\n" : "   ");
    }
  }

  if (os_ && i%5 != 0)
    *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doFinishOutput()
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
// Function      : FrequencyProbe::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyProbe::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doSteppingComplete()
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
