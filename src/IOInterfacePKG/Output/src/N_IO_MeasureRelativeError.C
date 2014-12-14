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
// Filename      : $RCSfile: N_IO_MeasureRelativeError.C,v $
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.17.2.1 $
// Revision Date  : $Date: 2014/08/28 22:41:46 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureRelativeError.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>
#include <N_DEV_Interpolators.h>
#include <Epetra_SerialDenseVector.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : RelativeError::RelativeError()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
RelativeError::RelativeError( const Util::OptionBlock & measureBlock):
  Base(measureBlock)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;
  if( comparisonFunctionName_.empty() )
    comparisonFunctionName_ = "L2NORM";

  if( !dataFileName_.empty() ) 
  {
    // load the data file
    std::ifstream dataFileStream( dataFileName_.c_str(), std::ios::in );
    if (!dataFileStream.good())
    {
      Report::UserError0() << "Could not open file in the measure function \"" + name_ + "\". Failed filename was = \"" + dataFileName_ + "\"";
    }   

    bool result = readData( dataFileStream, varNames_, dataValues_ );
    if( !result ) 
    {
      Report::UserError0() << "An error was encountered when reading the file in the measure function \"" << name_ << "\". Failed filename was = \"" << dataFileName_ << "\"";
    }
    dataFileStream.close();

    // copy over data into column vector for independent var
    indepVarValues_.resize( dataValues_.size() );
    double minIndepVarValues = 0.0;
    int maxVarColumn = ( independentVarColumn_ > dependentVarColumn_ ) ? independentVarColumn_ : dependentVarColumn_;
    for ( unsigned int i=0; i< dataValues_.size(); i++ )
    {
      int rowLength = dataValues_[i].size();

      if( maxVarColumn >= rowLength )
      {
        Report::UserError0() << "In measure function \"" << name_ 
                             << "\", using data from filename \"" << dataFileName_ 
                             << "\". Requested column for independent variable " 
	                     << independentVarColumn_ << " or dependent variable " 
                             << dependentVarColumn_ << " does not exist in the data file for entry " 
                             << i << std::endl;
      }
      indepVarValues_[i] = dataValues_[i][independentVarColumn_];

      // Check that the data is strictly monotonically increasing, otherwise throw an error
      if ( i )
      {
        if ( indepVarValues_[i] < minIndepVarValues )
        {
          Report::UserError0() << "In measure function \"" << name_ 
                               << "\", using data from filename \"" << dataFileName_ 
                               << "\". Independent variables are not sorted in monotonically increasing order." 
                               << std::endl;
        }
        minIndepVarValues = indepVarValues_[i];
      }
      else
      {
        // Initialize minIndepVarValues to the first entry.
        minIndepVarValues = indepVarValues_[0];
      }
    } 
 
    // results from the simulation to compare to a column in dataValues; 
    simulationDataVals_.reserve( dataValues_.size() );
    timeFreqVals_.reserve( dataValues_.size() );
  }
}

void RelativeError::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();
  
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for relative error measure, \"" + name_ + "\" Exiting.";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg);
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : RelativeError::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void RelativeError::reset() 
{
  initialized_=false;
  calculationDone_=false;
  simulationDataVals_.clear();
}


//-----------------------------------------------------------------------------
// Function      : RelativeError::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RelativeError::updateTran(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{  
  if( !calculationDone_ && withinFromToWindow( circuitTime ) )
  {
    timeFreqVals_.push_back( circuitTime );

    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, stateVec, storeVec, 0);
      simulationDataVals_.push_back( outVarValues_[i] );
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : RelativeError::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void RelativeError::updateDC(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec)
{
  //Xyce::dout() << "in RelativeError updateDC " << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : RelativeError::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2014
//-----------------------------------------------------------------------------
void RelativeError::updateAC(Parallel::Machine comm, const double frequency, const N_LAS_Vector * solnVec, const N_LAS_Vector *imaginaryVec)
{
  if( !calculationDone_ && withinFromToWindow( frequency ) )
  {
    timeFreqVals_.push_back( frequency );

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0, imaginaryVec);
      simulationDataVals_.push_back( outVarValues_[i] );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : RelativeError::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double RelativeError::getMeasureResult()
{
  // make an epetra vector to hold difference values
  Epetra_SerialDenseVector differenceVector(dataValues_.size());

  // Interpolate simulation data to input data values.
  std::vector<double> interpDataVals( dataValues_.size(), 0.0 );
  if (timeFreqVals_.size() > 0)
  {
    Xyce::Device::akima<double> interp;
    interp.init( timeFreqVals_, simulationDataVals_ );
    for (unsigned int i=0; i < dataValues_.size(); i++)
    {
      interp.eval( timeFreqVals_, simulationDataVals_, indepVarValues_[i], interpDataVals[i] );
    }
  }
  //Xyce::dout() << "Found " << numberPointsFound << " out of a possible " << dataValues.size() << std::endl;
  
  // load the difference vector
  for (unsigned int i=0 ; i<dataValues_.size(); i++)
  {
    differenceVector[i] = (interpDataVals[i]-dataValues_[i][dependentVarColumn_]);
  }

  if( comparisonFunctionName_ == "L1NORM" )
  {
    calculationResult_ = differenceVector.Norm1();
  }
  else if( comparisonFunctionName_ == "L2NORM" )
  {
    calculationResult_ = differenceVector.Norm2();
  }
  else
  {
    calculationResult_ = differenceVector.NormInf();
  }

  return calculationResult_;
}


//
// utility functions for reading in a file 
//

bool readData( std::istream & inputStream, 
               std::vector<std::string> & varNames, 
               std::vector< std::vector< double > > & dataValues )
{
  // read in the first line and try to devine the format of the file.
  // some possibliities are:
  //
  // (1) white space deliminated table with the first line optionally being 
  //     the name of each column and potentially a line of text 
  //     at the end.  This is the same as Xyce's prn format.
  //     Example: 
  //
  //     [ text  text  text ]
  //     number number number 
  //     ... 
  //     number number number 
  //     [ optional end of simulation text ] 
  //
  // (2) a compressed DC sweep format where the first line is the 
  //     name of the first independent variable and the second line
  //     gives the name of the second independent variable and 
  //     the values of the seconde indepentent variable for each column 
  //
  //     Example:
  //     text
  //     text value1 ... valueN
  //     value0 ... valueN
  //     value0 ... valueN
  //     value0 ... valueN
  //
  //     real example
  //     VD
  //     VG   -1.5 -0.5 0.0 0.5 1.0
  //     -0.9 0.00  0.1 0.2 0.1 0.1
  //     -0.8 0.00  0.1 0.2 0.1 0.1
  //
  // (3) probe or "CSD" formatted.  Starts with "#H" for header and "#N" for data


  // current line being parsed
  std::string workingLine;

  // try to read the first line
  std::getline( inputStream, workingLine );

  // a stringstream used in parsing up the line that was just read
  std::stringstream streamForWorkingLine;
  // store up data read for each line so it can 
  // later be saved in the vector< vector< double > >
  std::vector<double> aDataLine;

  
  while( inputStream.good() )
  { 
    // if fist line is "#H" then we know we have a probe formatted file
    if( workingLine.find( "#H" ) == 0 )
    {
      // have a probe formatted file so read that 
    }
  
    // look for one or more text headers (there may be none)
    // trying to decern "text" from "number" at this stage 
    //streamForWorkingLine.flush();
    streamForWorkingLine.clear();
    streamForWorkingLine << workingLine;
    std::string  subString;
    // break up line by whitespace 
    while( streamForWorkingLine >> std::ws >> subString )
    {
      if( Util::isValue( subString ) )
      {
        // found a numeric value 
        aDataLine.push_back( Util::Value(subString) );
      }
      else 
      {
        varNames.push_back(subString);
      }
    }
    
    if( !aDataLine.empty() )
    {
      dataValues.push_back( aDataLine );
      aDataLine.clear();
    }
    // 
    std::getline( inputStream, workingLine );
  } 

/*
  // for debugging purposes output the data read 
  Xyce::dout() << "Data read from stream " << std::endl;
  for( int i=0; i<varNames.size(); i++ )
    Xyce::dout() << " varNames[" << i << "]=" << varNames[i] << std::endl;

  Xyce::dout() << "Data values: " << std::endl;
  for( int i=0; i<dataValues.size(); i++ )
  {
    for( int j=0; j<dataValues[i].size(); j++ )
      Xyce::dout() << dataValues[i][j] << "\t";
    Xyce::dout() << std::endl;
  }
*/

  return true;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce
