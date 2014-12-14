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
// Filename       : $RCSfile: N_IO_OutputterFrequencyRaw.C,v $
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
// Revision Number: $Revision: 1.5.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterFrequencyRaw.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : FrequencyRaw
// Purpose       : Outputter class for frequency-domain output, rawfile output
//                 format
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::FrequencyRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRaw::FrequencyRaw(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    os_(0),
    numPoints_(0),
    numPointsPos_(0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::~FrequencyRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRaw::~FrequencyRaw()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyRaw::frequencyHeader(Parallel::Machine comm)
{
  if (os_)
  {
    std::ostream &os = *os_;

    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      os << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      os << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin();
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      plotName << "Transient Analysis";
    }
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {
      plotName << "AC Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    os << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    os << "Flags: " << flags << std::endl;

    int numVars = 0;
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {}
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {}
    else
    {
      ++numVars;
    }
    numVars += opList_.size();

    // format number of internal and external variables included here + time
    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = os_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the
      // in the raw file header.  Optionally let one output the version
      // if whatever program is going to read the file expects it
      os << "Version: " << Util::Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    // add timestep header info(Spice3f5 style)
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
    }
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {
    }
    else
    {
      os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
      ++i;
    }

    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      std::string tmpNodeName, tmpType;
      // set the type
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if ((*it)->getName() == "INDEX")  { }
      else if ((*it)->getName() == "TIME")  { tmpType = "time"; }
      else if ((*it)->getName() == "FREQUENCY")  { tmpType = "frequency"; }
      else if ((*it)->getName()[0] == 'I')  { tmpType = "current";    }
      else if ((*it)->getName()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*it)->getName()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyRaw::doOutputFrequency(
  Parallel::Machine     comm,
  double                frequency,
  const N_LAS_Vector &  real_solution_vector,
  const N_LAS_Vector &  imaginary_solution_vector)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, printParameters_.defaultExtension_, printParameters_.suffix_, outputManager_.getNetListFilename());
    os_ = outputManager_.openBinaryFile(outFilename_);

    numPoints_ = 0;
  }

  if (numPoints_ == 0)
    frequencyHeader(comm);

  // select values to write from .PRINT line if FORMAT=RAW
  for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    complex result = getValue(comm, *(*it), Util::Op::OpData(0, &real_solution_vector, &imaginary_solution_vector, 0, 0, 0));
    if (os_)
    {
      double realPart=result.real();
      double imagPart=result.imag();
      os_->write((char *) &realPart, sizeof( double));
      os_->write((char *) &imagPart, sizeof( double));
    }
  }

  // keep track of number of datapoints
  ++numPoints_;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyRaw::doFinishOutput()
{
  if (os_)
  {
    if (numPoints_ != 0)
    {
      // need to move file pointer back to header and
      // write out the number of points.
      long currentFelePost = os_->tellp();

      // locate the position for number of points
      os_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*os_) << numPoints_;

      // move file pointer to the end again.
      os_->seekp( currentFelePost);
    }
  }

  // reset numPoints_ as it is used as a flag to print the header.
  numPoints_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
