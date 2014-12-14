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
// Filename       : $RCSfile: N_IO_OutputMacroResults.C,v $
//
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.1 $
//
// Revision Date  : $Date: 2014/08/28 22:41:47 $
//
// Current Owner  : $Author: rlschie $
//-------------------------------------------------------------------------

#include <sstream>
#include <Xyce_config.h>

#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

#include <N_IO_Objective.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_FourierMgr.h>

#include <Teuchos_oblackholestream.hpp>

namespace Xyce {
namespace IO {

typedef std::map<std::string, Objective> ObjectiveMap;
typedef std::vector<std::pair<std::string, std::string> > StringPairVector;

//-----------------------------------------------------------------------------
// Function      : outputMacroResults
// Purpose       : if any post-process analysis was specified(like objective or measure)
//                 output results.  Called after simulation is over.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL Electrical and Microsystem Modeling
// Creation Date : 11/05/08
//-----------------------------------------------------------------------------
void outputMacroResults(
  Parallel::Machine             comm,
  const ObjectiveMap &          objective_map,
  const Measure::Manager &      measure_manager,
  FourierMgr &                  fourier_manager,
  std::string                   netlist_filename,
  const StringPairVector &      response_functions,
  std::string                   outputResponseFilename,
  const int                     step_number)
{
  if (!objective_map.empty())
  {
    Xyce::dout() << std::endl
      << " ***** Analysis Functions ***** " << std::endl
      << std::endl;

    for (ObjectiveMap::const_iterator it = objective_map.begin(), end = objective_map.end(); it != end ; ++it)
    {
      double val = (*it).second.evaluate();
      Xyce::dout() << (*it).first << " = " << val << std::endl;
    }
  }

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results only if .measure is being performed on any variables.
  if (measure_manager.isMeasureActive())
  {
    // Output the Measure results to Xyce::dout().
    // Make sure the function gets called on all processors, but only one outputs it.
    if (Parallel::rank(comm) == 0)
    {
      measure_manager.outputResults( Xyce::dout() );
    }
    else
    {
      measure_manager.outputResults( outputBHS );
    }

    // Output the Measure results to file.
    if (Parallel::rank(comm) == 0)
    {
      // this is a hack until this is function is refactored.
      // If this is a ".step" analysis then the latest measure
      // results have already been written.  So the following code 
      // with just overwrite the file with the same results.  
      if( step_number == 0 )
      {
        std::ostringstream converterBuff;
        converterBuff << step_number;
        std::string filename = netlist_filename + ".mt" + converterBuff.str();
        outputFileStream.open( filename.c_str() );
        measure_manager.outputResults( outputFileStream, false );
        outputFileStream.close();
      }
    }
    else {
      measure_manager.outputResults( outputBHS, false );
    }
  }

  // Output the Fourier results to file if Fourier analysis is being performed.
  // Make sure the function gets called on all processors, but only one outputs it.
  if (fourier_manager.isFourierActive())
  {
    if (Parallel::rank(comm) == 0)
    {
      std::string filename = netlist_filename + ".four";
      outputFileStream.open( filename.c_str() );
      fourier_manager.outputResults(outputFileStream);
      outputFileStream.close();
    }
    else {
      fourier_manager.outputResults( outputBHS );
    }
  }

  // if the response list is not empty, try to dump those results to a file
  // a big limitation here is that all responses must be measure functions.
  // need to make this more flexible to include all solution vars and
  // objectives and results.
  if (!response_functions.empty())
  {
    std::ofstream responseOFS(outputResponseFilename.c_str());

    for (StringPairVector::const_iterator currRespItr = response_functions.begin(), endRespItr = response_functions.end();
         currRespItr != endRespItr; ++currRespItr)
    {
      double respvalue = -1.0;
      // need to parse name from XXX_X:name So find last ":" and extract
      // remaining string as value that dakota is looking for.
      std::string tempName = currRespItr->first;
      std::string::size_type beginningOfVarName = tempName.find_last_of(":");

      if (beginningOfVarName != std::string::npos)
      {
        int numChars = (currRespItr->first).length() - beginningOfVarName;
        tempName.assign(currRespItr->first, beginningOfVarName + 1, numChars);
      }
      //Xyce::dout() << " looking for " << currRespItr->first << std::endl;
      ExtendedString es(tempName);
      bool found = false;
      measure_manager.getMeasureValue(es.toUpper(), respvalue, found);
      if (!found)
      {
        // look for value in .objective functions
        for (ObjectiveMap::const_iterator it = objective_map.begin(), end = objective_map.end(); it != end ; ++it)
        {
          ExtendedString name((*it).first);
          if( es.toUpper() == name.toUpper() )
          {
            respvalue = (*it).second.evaluate();
            found = true;
          }
        }
      }
      responseOFS << respvalue   << "   " << tempName << std::endl;
    }
    responseOFS.close();
  }
}

} // namespace IO
} // namespace Xyce
