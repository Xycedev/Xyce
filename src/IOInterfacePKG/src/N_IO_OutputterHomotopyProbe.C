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
// Filename       : $RCSfile: N_IO_OutputterHomotopyProbe.C,v $
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
// Revision Number: $Revision: 1.6.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHomotopyProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : HomotopyProbe
// Purpose       : Outputter class for homotopy output, Probe output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyProbe::HomotopyProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeHomotopy_(true)
{
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::~HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyProbe::~HomotopyProbe()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::homotopyHeader(
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const N_LAS_Vector &                  solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  index_ = 0;

  printCount_ = opList_.size();

  os << "#H" << std::endl;

  os << "SOURCE='Xyce' VERSION='"
                   << Util::Version::getShortVersionString() << "'" << std::endl;

  os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;
  os << "SUBTITLE='spice probe data'" << std::endl;

  // set the time/date stamp
  os << getTimeDateStamp();
  os.setf(std::ios::scientific);
  os.precision(printParameters_.streamPrecision_);
  os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
  os << "'" << std::endl;

  if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    os << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
  else
    os << "ANALYSIS='DC transfer characteristic' " << "SERIALNO='12345'" <<  std::endl;

  os << "ALLVALUES='NO' COMPLEXVALUES='NO' " << "NODES='" << printCount_ << "'" << std::endl;

  if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
  {
    os << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" << std::endl;
  }
  else
  {
    os << "SWEEPVAR='Voltage' SWEEPMODE='VAR_STEP'" <<std::endl;
  }

  // This line assumes that we're doing a homotopy that goes from 0 to 1.
  // This will never be a transient output.
  os << "XBEGIN='0.0'  XEND='1.0'"<<std::endl;

  os << "FORMAT='0 VOLTSorAMPS;EFLOAT : " << "NODEorBRANCH;NODE  '  " << std::endl;

  os << "DGTLDATA='NO'" << std::endl;

  os << "#N" << std::endl;

  // print the continuation parameter names:
  int i = 0;
  for (std::vector<std::string>::const_iterator iter_name = parameter_names.begin(); iter_name != parameter_names.end(); ++iter_name, ++i)
  {
    os << "'" << *iter_name << "' ";

    if (i > 3)
    {
      i = 0;
      os << std::endl;
    }
  }

  // print output variable names:
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++i)
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


  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const N_LAS_Vector &                  solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  double tmpTime = outputManager_.getCircuitTime(); // outputManager_.getAnaIntPtr()->getTime();

  if (Parallel::rank(comm) == 0)
  {

    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());

      if (Parallel::rank(comm) == 0 && outStreamPtr_ == 0)
      {
        outStreamPtr_ = outputManager_.openFile(outFilename_);
      }
      homotopyHeader(parameter_names, parameter_values, solution_vector);

      firstTimeHomotopy_ = false;
    }

    os.width( 0);
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      os << "#C " << tmpTime << " ";
      os << printCount_ << std::endl;
    }
    else
    {
      os << "#C " << outputManager_.getPRINTDCvalue() << " ";
      os << printCount_ << std::endl;
    }

    //-------------------------------------
    //HOMOTOPY PARAM VALUE OUTPUT GOES HERE
    //-------------------------------------

    for (int iparam=0;iparam < parameter_values.size(); ++iparam)
    {

      if (printParameters_.delimiter_ == "")
      {
        os.width(printParameters_.streamWidth_);
      }
      else
      {
        os.width(0);
        if (printParameters_.delimiter_ != "")
          os << printParameters_.delimiter_;
      }

      os << parameter_values[iparam];
    }
  } // procID

  Util::Op::OpList::const_iterator iterParam = opList_.begin();
  Util::Op::OpList::const_iterator last = opList_.end();

  int i;
  for (i = 1; iterParam != last; ++iterParam, ++i)
  {
    double result = getValue(comm, *(*iterParam), Util::Op::OpData(0, &solution_vector, 0, 0, 0, 0)).real();
    if (Parallel::rank(comm) == 0)
    {
      if (printParameters_.delimiter_ == "")
      {
        os.width(printParameters_.streamWidth_);
      }
      else
      {
        os.width(0);
        if (printParameters_.delimiter_ != "")
          os << printParameters_.delimiter_;
      }
      os << result;
    }
  }

  if (Parallel::rank(comm) == 0)
    os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyProbe::doFinishOutput()
{
  firstTimeHomotopy_ = true;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doStartStep(
  int                           step,
  int                           max_step)
{
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doSteppingComplete()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    (*outStreamPtr_) << "#;" << std::endl;
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
