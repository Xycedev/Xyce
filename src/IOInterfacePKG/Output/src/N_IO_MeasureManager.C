//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_IO_MeasureManager.C,v $
// Purpose       : This file contains the functions to manage measure objects
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.34.2.4 $
// Revision Date  : $Date: 2014/09/05 17:24:59 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <utility>
#include <sstream>

#include <Teuchos_oblackholestream.hpp>

#include <N_IO_OutputPrn.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_MeasureRiseFallDelay.h>
#include <N_IO_MeasureAverage.h>
#include <N_IO_MeasureMax.h>
#include <N_IO_MeasureMin.h>
#include <N_IO_MeasurePeakToPeak.h>
#include <N_IO_MeasureRMS.h>
#include <N_IO_MeasureFrequency.h>
#include <N_IO_MeasureDuty.h>
#include <N_IO_MeasureOnTime.h>
#include <N_IO_MeasureOffTime.h>
#include <N_IO_MeasureFindWhen.h>
#include <N_IO_MeasureEquationEvaluation.h>
#include <N_IO_MeasureDerivativeEvaluation.h>
#include <N_IO_MeasureIntegralEvaluation.h>
#include <N_IO_MeasureRelativeError.h>
#include <N_IO_MeasureFourier.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_PDS_ParMap.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Class         : MeasureOptionsReg
// Purpose       : functor for registering Measure options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct MeasureOptionsReg : public PkgOptionsReg
{
  MeasureOptionsReg( Measure::Manager &mgr )
    : measureManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return measureManager_.addMeasure( options ); }

  Measure::Manager &measureManager_;
};

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool Manager::registerPkgOptionsMgr(PkgOptionsMgr &pkgOpt)
{
  pkgOpt.submitRegistration(
    "MEASURE", netListFilename_, new MeasureOptionsReg(*this));
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::Manager(
  const std::string &   netlist_filename)
  : netListFilename_(netlist_filename)
{}

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::~Manager()
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    delete (*it);
}

void
Manager::makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager) 
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    (*it)->makeMeasureOps(comm, op_builder_manager);
}

void
Manager::notify(
  const Analysis::StepEvent &   step_event)
{
  switch (step_event.state_) {
    case Analysis::StepEvent::INITIALIZE:
      for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
      {
        (*it)->reset();
      }
      break;

    case Analysis::StepEvent::STEP_STARTED:
      for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
      {
        (*it)->reset();
      }
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      {
        if( !allMeasuresList_.empty() )
        {
          std::ostringstream converterBuff;
          converterBuff << step_event.count_;
          std::string filename = netListFilename_ + ".mt" + converterBuff.str();
          std::ofstream outputFileStream;
          outputFileStream.open( filename.c_str() );
          outputResults( outputFileStream, false );
          outputFileStream.close();
        }
      }

      break;

    case Analysis::StepEvent::FINISH:
      break;
  }
}


//-----------------------------------------------------------------------------
// Function      : Manager::addMeasure
// Purpose       : Entery point when .measure lines are pased in the netlist
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Manager::addMeasure(const Util::OptionBlock & measureBlock )
{
  // based on what's in the option block passed in, we
  // create the needed measure instance
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "In Manager::addMeasure" << std::endl;
  Xyce::dout() << "measureLine passed it was: " << std::endl << measureBlock << std::endl;
#endif

  // here's the base data we should pull from the option block while
  std::string type;

  if( measureBlock.tagExists( "TYPE" ) )
  {
    type = measureBlock.getTagValueAsString("TYPE");
  }
  else
  {
    // this shouldn't happen, but catch it if does
    Report::UserError0() << "Missing TYPE in Manager";
  }

  N_IO_MeasureBase  * theMeasureObject = 0;
  // have enough info to make the correct measure class
  if( type=="TRIG" || type=="TARG" )
  {
    theMeasureObject = new N_IO_MeasureRiseFallDelay( measureBlock);
  }
  else if( type=="AVG" )
  {
    theMeasureObject = new N_IO_MeasureAverage( measureBlock);
  }
  else if( type=="MAX" )
  {
    theMeasureObject = new N_IO_MeasureMax( measureBlock);
  }
  else if( type=="MIN" )
  {
    theMeasureObject = new N_IO_MeasureMin( measureBlock);
  }
  else if( type=="PP" )
  {
    theMeasureObject = new N_IO_MeasurePeakToPeak( measureBlock);
  }
  else if( type=="RMS" )
  {
    theMeasureObject = new N_IO_MeasureRMS( measureBlock);
  }
  else if( type=="FREQ" )
  {
    theMeasureObject = new N_IO_MeasureFrequency( measureBlock);
  }
  else if( type=="DUTY" )
  {
    theMeasureObject = new N_IO_MeasureDuty( measureBlock);
  }
  else if( type=="ON_TIME" )
  {
    theMeasureObject = new N_IO_MeasureOnTime( measureBlock);
  }
  else if( type=="OFF_TIME" )
  {
    theMeasureObject = new N_IO_MeasureOffTime( measureBlock);
  }
  else if( type=="FIND" || type=="WHEN" )
  {
    theMeasureObject = new N_IO_MeasureFindWhen( measureBlock);
  }
  else if( type=="PARAM" || type=="EQN"  )
  {
    theMeasureObject = new N_IO_MeasureEquationEvaluation( measureBlock);
  }
  else if( type=="DERIVATIVE" || type=="DERIV" )
  {
    theMeasureObject = new N_IO_MeasureDerivativeEvaluation( measureBlock);
  }
  else if( type=="INTEGRAL" || type=="INTEG" )
  {
    theMeasureObject = new N_IO_MeasureIntegralEvaluation( measureBlock);
  }
  else if( type=="ERROR" )
  {
    theMeasureObject = new N_IO_MeasureRelativeError( measureBlock);
  }
  else if( type=="FOUR" )
  {
    theMeasureObject = new N_IO_MeasureFourier( measureBlock);
  }
  else
  {
    // unknown type issue warning.
    Xyce::Report::UserWarning() << "Unknown MEASURE type requested, \"" << type << "\".  This measure will be ignored";
  }

  // if the measure object is supported, then add it to the active and all lists
  if (theMeasureObject && theMeasureObject->typeSupported_ )
  {
    allMeasuresList_.push_back( theMeasureObject );
    activeMeasuresList_.push_back( theMeasureObject );
  }
  else
    delete theMeasureObject;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateTranMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateTranMeasures(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  // loop over active masure objects and get them to update themselves.
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateTran(comm, circuitTime, solnVec, stateVec, storeVec);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end());
}


//-----------------------------------------------------------------------------
// Function      : Manager::updateDcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateDCMeasures(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateDC(comm, dcParamsVec, solnVec, stateVec, storeVec);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end());
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateAcMeasures
// Purpose       : Called during the simulation to update the measure objects
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/21/2014
//-----------------------------------------------------------------------------
void Manager::updateACMeasures(Parallel::Machine comm, const double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateAC(comm, frequency, real_solution_vector, imaginary_solution_vector);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), std::mem_fun(&Measure::Base::finishedCalculation)),
                            activeMeasuresList_.end()); }

//-----------------------------------------------------------------------------
// Function      : Manager::outputResults
// Purpose       : Output measure results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
std::ostream &Manager::outputResults( std::ostream& outputStream, bool printHeader ) const
{
  if (!allMeasuresList_.empty())
  {
    if (printHeader)
    {
      outputStream << std::endl
                   << " ***** Measure Functions ***** " << std::endl
                   << std::endl;
    }

    // loop over active measure objects and get the results
    for (MeasurementVector::const_iterator it = allMeasuresList_.begin(), end = allMeasuresList_.end(); it != end; ++it)
    {
      (*it)->printMeasureResult( outputStream );
    }
  }
  
  return outputStream;
}


//-----------------------------------------------------------------------------
// Function      : Manager::getMeasureValue
// Purpose       : Get the value of a .measure test
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2010
//-----------------------------------------------------------------------------
void Manager::getMeasureValue(const std::string &name, double &value, bool &found) const
{
  found = false;

  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
  {
    if (equal_nocase((*it)->name_, name))
    {
      value = (*it)->getMeasureResult();
      found = true;
      return;
    }
  }
}

const Base *Manager::find(const std::string &name) const
{
  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
  {
    if (equal_nocase((*it)->name_, name))
    {
      return *it;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::remeasure
// Purpose       : Recompute measure functions based on existing Xyce output
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/8/13
//-----------------------------------------------------------------------------
void Manager::remeasure(N_PDS_Comm &pds_comm, const std::string &netlist_filename, const std::string &remeasure_path, bool transient_flag, OutputMgr &output_manager)
{
  Xyce::dout() << "In OutputMgr::remeasure " << std::endl
               << "file to reprocess through measure functions. " << remeasure_path << std::endl;

  // open file for reading
  // just support PRN format for now
  OutputFileBase * fileToReMeasure = new N_IO_OutputPrn();
  if (!fileToReMeasure->openFileForRead(remeasure_path))
  {
    // open failed.  report error and exit remeasure routine

    return;
  }

  // load data-names & allocate space for a line
  std::vector<std::string> fileVarNames;
  if (!(fileToReMeasure->getOutputVarNames(fileVarNames)))
  {
    // reading var names failed.  report error and exit remeasure routine.

  }

  // this code will need to move to the output handler classes as it will be
  // specific to the file type(i.e. some formats store the entire solution
  // regardless of what is on the print line.
  fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames);

//   Xyce::dout() << "Original var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     Xyce::dout() << "\"" << fileVarNames[i] << "\", ";
//   }
//   Xyce::dout() << std::endl;
//
//   fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames);
//
//   Xyce::dout() << "extracted var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     Xyce::dout() << "\"" << fileVarNames[i] << "\", ";
//   }
//   Xyce::dout() << std::endl;

  // assume analysis type is DC(Xyce hasn't processed the analysis type yet
  // when this function is called.  In the next loop if we find a column of
  // data for TIME then we can treat this as a transient analysis.

  int timeIndex=0;
  // set up allNodes map to map var names to indices in the solution vector
  int numVars = fileVarNames.size();
  for (int i=0; i<numVars; i++)
  {
    output_manager.getAllNodes()[fileVarNames[i]] = std::make_pair(i, 0);
    ExtendedString tmpStr(fileVarNames[i]);
    tmpStr.toUpper();
    output_manager.getAllNodes()[tmpStr] = std::make_pair(i, 0);
    // while scanning fileVarNames look for "TIME" as a name in the 0th or 1st
    // column.  We will use this as a key to figure out if the output file is
    // transient or DC data(no support for .measure in other analysis types yet)
    if ((i<2) &&(fileVarNames[i]=="TIME"))
    {
      timeIndex = i;
    }
  }

  // create an N_LAS_Vector to hold the data from the file.  This is
  // needed to support getPrintValue_ for evaluating the measure functions.
  std::vector<int> lbMap;
  lbMap.resize(numVars);
  for (int i=0;i<numVars; i++)
  {
    lbMap[i]=i;
  }

  {
    N_PDS_ParMap aParMap(numVars, numVars, lbMap, 0, &pds_comm);
    N_LAS_Vector varValuesVec(aParMap);

    varValuesVec.putScalar(0);

    // set up context for items in .measure line

    // run though lines in the file calling update measure as we go.
    while (fileToReMeasure->getOutputNextVarValues(&varValuesVec))
    {
      if (transient_flag)
      {
        output_manager.setCircuitTime(varValuesVec[timeIndex]);
        updateTranMeasures(pds_comm.comm(), output_manager.getCircuitTime(), &varValuesVec, 0, 0 );
      }
      else
      {
        updateDCMeasures(pds_comm.comm(), output_manager.getDCParamVec(), &varValuesVec, 0, 0 );
      }
      varValuesVec.putScalar(0);
    }
  }

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;

  // Output the Measure results to Xyce::dout().
  // Make sure the function gets called on all processors, but only one outputs it.
  if (Parallel::rank(pds_comm.comm()) == 0)
  {
    outputResults( Xyce::dout() );
  }
  else
  {
    outputResults( outputBHS );
  }

  // Output the Measure results to file.
  if (Parallel::rank(pds_comm.comm()) == 0)
  {
    // Adding "0" to the end of this string for the step number, which was always 0
    // for .measure previously anyways.
    std::string filename = netlist_filename + ".mt0";
    std::ofstream outputFileStream(filename.c_str());
    outputResults(outputFileStream );
    outputFileStream.close();
  }
  else {
    outputResults( outputBHS );
  }

  fileToReMeasure->closeFileForRead();
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
