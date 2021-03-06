//-----------------------------------------------------------------------------
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
// Filename       : $N_DAK_DakotaInterface.C$
//
// Purpose        : This class defines an interface between the optimization
//                  and analysis routines in Dakota and the circuit routines
//                  in Xyce.  The object here is allow Xyce to read a netlist
//                  with some Dakota commands embeded within the netlist and
//                  pass off control to Dakota to organize and run as many
//                  Xyce simulations as are needed for Dakota's analysis.
//
// Special Notes  :
//
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.40 $
//
// Revision Date  : $Date: 2014/07/15 19:03:33 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DAK_DakotaInterface.h>

// Xyce includes
#include <N_CIR_Xyce.h>
#include <N_DEV_DeviceInterface.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_System.h>
#include <N_ANP_AnalysisManager.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Misc.h>

#include <sstream>

#include <iostream>

N_DAK_DakotaInterface::N_DAK_DakotaInterface( const Dakota::ProblemDescDB& problem_db):
  Dakota::DirectApplicInterface( problem_db )
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::N_DAK_DakotaInterface() " << std::endl;
#endif
}

N_DAK_DakotaInterface::~N_DAK_DakotaInterface()
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::~N_DAK_DakotaInterface() " << std::endl;
#endif
  deleteCargs( xyceIargs, xyceCargs );
}

//
// 
//
void N_DAK_DakotaInterface::setArguments( int iargsIn, char * cargsIn [] )
{
  // copy the iargs and cargs
  xyceIargs = iargsIn;
  copyCargs( xyceIargs, cargsIn, xyceCargs );
  
#ifdef Xyce_Dakota_Debug
	Xyce::dout() << "N_DAK_DakotaInterface::setArguments" << std::endl;
	for(int i=0; i<xyceIargs; i++)
	{
	  Xyce::dout() << "arg[" << i << "] = \"" << xyceCargs[i] << "\"" << std::endl;
	}
#endif
  
}

//
// routines called directly by Dakota
//

int N_DAK_DakotaInterface::derived_map_ac( const Dakota::String& ac_name)
{

#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::derived_map_ac() ac_name = \"" << ac_name << "\"" << std::endl;
  Xyce::dout() << "numACV = " << numACV << " numADIV = " << numADIV << " numADRV = " << numADRV << std::endl;
  Xyce::dout() << "xCLabels.size() = " << xCLabels.size() 
               << "  xDILabels.size() = " << xDILabels.size()
               << "  xDRLabels.size() = " << xDRLabels.size() << std::endl;
#endif
  //
  // loop over all the continuous and discrete variables labels and values 
  // 
  std::ostringstream varValue;
  // make sure the list of variable substitutions is clear
  variableSubVec.clear();
  
  // now fill the variable substitution list with fresh values
  // pack in the real values
  varValue.precision(16);
  for(int i=0; i < numACV; i++)
  {
    // clear the string stream
    varValue.str( std::string() );
    varValue.clear();
    varValue << xC[i];
    std::string varValueAsString( varValue.str() );
    variableSubVec.push_back( std::make_pair<std::string,std::string>( xCLabels[i], varValueAsString ) );
  }
	
  // now get discrete integer variables 
  varValue.precision(0);
  for(int i=0; i < numADIV; i++)
  {
    varValue.str( std::string() );
    varValue.clear();
    varValue << xDI[i];
    std::string varValueAsString( varValue.str() );
    variableSubVec.push_back( std::make_pair<std::string,std::string>( xDILabels[i], varValueAsString ) );
  }
  
  // now get discrete real variables 
  varValue.precision(0);
  for(int i=0; i < numADRV; i++)
  {
    varValue.str( std::string() );
    varValue.clear();
    varValue << xDR[i];
    std::string varValueAsString( varValue.str() );
    variableSubVec.push_back( std::make_pair<std::string,std::string>( xDRLabels[i], varValueAsString ) );
  }

#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::derived_map_ac() variableSubVec" << std::endl;
  std::vector< std::pair<std::string,std::string> >::iterator currPair = variableSubVec.begin();
  std::vector< std::pair<std::string,std::string> >::iterator endPair = variableSubVec.end();
  while( currPair != endPair )
  {
    Xyce::dout() << "\"" << currPair->first << "\", \"" << currPair->second << "\"" << std::endl;
    currPair++;
  }
#endif

  // create the Xyce object and pass it the list of variables and arguments
  xyceCirPtr_ = rcp( new N_CIR_Xyce() );
  xyceCirPtr_->setNetlistParameters( variableSubVec );
  xyceCirPtr_->initialize( xyceIargs, xyceCargs );

  //  
  // look at the response functions requested
  //
	 
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::derived_map_ac() Requested response functions: " 
               << numFns << std::endl; 
#endif

  initResponse();
   
  //
  // run the circuit
  // 	
  // now run the simulation 
  // in the future, there should be a try{ }catch() block around this.
  bool simResult = xyceCirPtr_->runSimulation();
	
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::derived_map_ac() Xyce iteration done." << std::endl;
#endif
 
  if( numFns )
  {
    // 
    // collect results from response functions
    //
    xyceCirPtr_->obtainResponses( fnLabels, simulationDataValues );
		
    for( int i=0; i<numFns; ++i )
    {
      fnVals[i] = simulationDataValues.at( i );
      if( !simResult )
      {
        // if the simulation failed, override the result by returning a nan
        // so the optimizer isn't confused by a zero result when there wasn't one
        const char * taqp=NULL;
        fnVals[i] = nan( taqp );
      }
    }
  }
	
  // 
  // now clean up the Xyce object
  //
  xyceCirPtr_->finalize();	
  xyceCirPtr_ = Teuchos::null;
  
  return 0;
}

void N_DAK_DakotaInterface::initResponse()
{
  // Check that the Response variables have been registered with .measure.
  //
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "N_DAK_DakotaInterface::initResponse() numFns = " << numFns << std::endl;
  for( int i=0; i<numFns; i++ )
  {
    Xyce::dout() << "Response variable " << i << ": \"" << fnLabels[i] << "\"" << std::endl;
  }
#endif

  // Size the endpoint simulation array and check that these function labels are valid.
  simulationDataValues.resize( numFns );
  bool success = xyceCirPtr_->checkResponseVars( fnLabels );

  if (!success)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Cannot locate all response variables required to execute Xyce-Dakota analysis");
  }
}


void N_DAK_DakotaInterface::copyCargs( const int originalIargs, char ** const originalCargs, char ** & copyCargs )
{
  copyCargs = new char * [originalIargs];
  if( copyCargs == NULL )
  {
    // remaining arguments aren't enough to get at least a netlist.
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Couldn't allocate memory for arguments in N_DAK_DakotaInterface::copyCargs().");
  }
  
  for (int i=0; i<originalIargs;++i)
  {
    // this seems odd, but it was how N_IO_CmdParse did this as well.
    if (originalCargs[i] == NULL) 
    {
      copyCargs[i] = NULL;
      continue;
    }
    
    std::string tmpString(originalCargs[i]);
    int size = tmpString.size()+2;
    copyCargs[i] = new char[size];
    if( copyCargs[i] == NULL )
    {
      // remaingin arguments aren't enough to get at least a netlist.
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Couldn't allocate memory for arguments in N_DAK_DakotaInterface::copyCargs().");
    }
    
    for (int j=0;j<size;++j)
    {
      copyCargs[i][j] = 0;
    }
    sprintf(copyCargs[i], tmpString.c_str());
  }
}

void N_DAK_DakotaInterface::deleteCargs( const int len, char ** & cargs )
{
  for(int i=0; i<len; i++)
  {
    if( cargs[i] )
    {
      delete [] cargs[i];
      cargs[i]=NULL;
    }
  }
  if( cargs )
  {
    delete cargs;
    cargs=NULL;
  }
}

